#!/usr/bin/perl -w

use strict;
use POSIX;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev', '/gpfs/fs121/h/toung/dev';
use File::Basename;
use PullGenomeSequence;
use Bio::SeqIO;
use Bio::Seq;

umask 0007;

my $rddFile;
my $tag;
my $region; ### (optional) e.g. 1-100-10000 ==> part = 1, starts at 100, ends at 10000
my $blat;
my $index;
my $fasta_folder; ### (optional)
my $pslx_folder; ### (optional)
my $flanking_distance; ### e.g. 25,50,75
my $kmer = 12;
my $stepSize = 5;
my $minIdentity;
my $repMatch = 2253;

my $options = GetOptions(
	"rddFile=s" => \$rddFile,
	"tag=s" => \$tag,
	"region=s" => \$region,
	"blat=s" => \$blat,
	"index=s" => \$index,
	"fasta_folder=s" => \$fasta_folder,
	"pslx_folder=s" => \$pslx_folder,
	"flanking_distance=s" => \$flanking_distance,
	"kmer=i" => \$kmer,
	"stepSize=i" => \$stepSize,
	"minIdentity=i" => \$minIdentity,
	"repMatch=i" => \$repMatch
);

(print STDERR "rddFile:\t$rddFile\n") && -e $rddFile or die "[STDERR]: rddFile doesn't exist\n";
(print STDERR "tag:\t$tag\n") && defined $tag or die "[STDERR]: tag doesn't exist\n";
my ($regionPart,$regionStart,$regionEnd);
if (defined $region) {
	($regionPart,$regionStart,$regionEnd) = split('-',$region);
	print STDERR "regionPart:\t$regionPart\n";
	print STDERR "regionStart:\t$regionStart\n";
	print STDERR "regionEnd:\t$regionEnd\n";
}
(print STDERR "blat:\t$blat\n") && -e $blat or die "[STDERR]: cant find blat\n";
(print STDERR "index:\t$index\n") && -e $index or die "[STDERR]: undefined $index\n";
if (defined $fasta_folder) {
	(print STDERR "fasta_folder:\t$fasta_folder\n") && -d $fasta_folder or die "[STDERR]: fasta folder doesn't exist\n";
}
if (defined $pslx_folder) {
	(print STDERR "pslx_folder:\t$pslx_folder\n") && -d $pslx_folder or die "[STDERR]: pslx folder doesn't exist\n";
}
(print STDERR "kmer:\t$kmer\n") && defined $kmer or die "[STDERR]: undefined kmer\n";
(print STDERR "stepSize:\t$stepSize\n") && defined $stepSize or die "[STDERR]: undefined stepSize\n";
(print STDERR "minIdentity:\t$minIdentity\n") && defined $minIdentity or die "[STDERR]: undefined minIdentity\n";
(print STDERR "repMatch:\t$repMatch\n") && defined $repMatch or die "[STDERR]: undefined repMatch\n";

my @flanking_distance = split(',',$flanking_distance);

foreach my $fl_dist (@flanking_distance) { 
	(print STDERR "flanking_distance:\t$fl_dist\n") && $fl_dist =~ /[0-9]+/ or die "[STDERR]: flanking distance invalid\n"; 
}

### make fasta files with sequences for each flanking distance
my %output;
foreach my $fl_dist (@flanking_distance) {
	my $fasta_file;
	$fasta_file = $fasta_folder . "/" if defined $fasta_folder;
	$fasta_file .= $tag . "_";
	$fasta_file .= $regionPart . "_" . $regionStart . "_" . $regionEnd . "_" if defined $region;
	$fasta_file .= "flankDist_" . $fl_dist . ".fasta";
	print STDERR "fasta_file:\t$fasta_file\n";
	$output{$fl_dist}{'fasta_file'} = $fasta_file;
	$output{$fl_dist}{'fasta_obj'} = Bio::SeqIO->new(-file => ">$output{$fl_dist}{'fasta_file'}", -format => 'fasta');
}

open(FILE,$rddFile) or die "[STDERR]: can't open $rddFile: $!\n";

while(<FILE>) {
	chomp;

	next if (/^#/);
	
	my @split = split('\t');
	
	my $lineNumber = $.;

	if (defined $region) {
		next if ($. < $regionStart);
		next if ($. > $regionEnd);
	}
			
	make_sequences(\@split,\@flanking_distance,\%output,$index,$lineNumber);
}
close(FILE);

## RUN BLAT ON EACH OF THE FASTA FILES ########################################
my %blat_results;
foreach my $fl_dist (@flanking_distance) {
	my $pslx;
	$pslx = $pslx_folder . "/" if defined $pslx_folder;
	$pslx .= $tag . "_";
	$pslx .= $regionPart . "_" . $regionStart . "_" . $regionEnd . "_" if defined $region;
	$pslx .= "flankDist_" . $fl_dist . ".pslx";

	print STDERR "pslx_file:\t$pslx\n";

	my $blat_command = "time $blat -stepSize=$stepSize -minIdentity=$minIdentity -repMatch=$repMatch -noHead -out=pslx $index $output{$fl_dist}{'fasta_file'} $pslx";
	print STDERR $blat_command, "\n";
	!system($blat_command) or die "[STDERR]: can't run $blat_command: $!\n";
	print STDERR "completed blatting for fl_dist $fl_dist\n";
}

print STDERR "completed BlatRDDFlankingSeq\n";

#######################################

sub make_sequences {
	my ($split,$flank_dist,$output,$index,$lineNumber) = @_;
	
	my $chrom = $split->[0];
	my $start = $split->[1];
	my $end = $split->[2];
	my $strand = $split->[3];
	
	($strand eq '+' || $strand eq '-' || $strand eq '+,-') or die "[STDERR]: weird strand $strand\n";
	
	my $ref_base = $split->[4];
	my $rdd_base = $split->[5];
	
	if ($strand eq '-') {
		$ref_base =~ tr/ACGT/TGCA/;
		$rdd_base =~ tr/ACGT/TGCA/;
	}

	foreach my $fl_dist (@{$flank_dist}) {
		my $display_id = join(";",$chrom,$end,$strand,$fl_dist,$split->[4],$split->[5],$lineNumber);
		
		my $region = $chrom . ":" . ($end - $fl_dist) . "-" . ($end + $fl_dist);
		my $sequence = pull_sequence($region,$index);
		
		### sequence may not be defined if at ends of chromosome
		next unless $sequence && length($sequence) == ($fl_dist*2 + 1);
		
		my @sequence = split('',$sequence);
		my @RNA_sequence = @sequence;
		if ($rdd_base eq 'X') {
			$RNA_sequence[$fl_dist] = '';	
		} else {
			$RNA_sequence[$fl_dist] = $rdd_base;
		}
		my $RNA_sequence = join('',@RNA_sequence);
		my $rna_obj = Bio::Seq->new(-seq => $RNA_sequence, -display_id => $display_id, -desc => "", -alphabet => "dna");
		$output->{$fl_dist}{'fasta_obj'}->write_seq($rna_obj);
	}
}


