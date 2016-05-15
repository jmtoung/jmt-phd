#!/usr/bin/perl -w

use strict;
use POSIX;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev', '/gpfs/fs121/h/toung/oldhome/dev';
use File::Basename;
use PullGenomeSequence;
use Bio::SeqIO;
use Bio::Seq;

umask 0007;

my $rddFile = "/gpfs/fs121/h/toung/work/2012-04-24_WholeFlowCell/step1-PoolRDD/step5-FilterOnTrueCounts/gsnap.GM12004_GM12750.defv20120410.unique.wholeFlowCellpooledFinal.min5pct2rddReads10totalReads.rdd";
my $tag = "test";
my $region = "1-1-1000"; ### (optional) e.g. 1-100-10000 ==> part = 1, starts at 100, ends at 10000
my $bam = "/gpfs/fs121/h/toung/aln/baseline/GM12750_deeplane1_trimlowqual_trimadapterv2/gsnap.GM12750_deeplane1_trimlowqual_trimadapterv2.defv20120410.unique.bam,/gpfs/fs121/h/toung/aln/baseline/GM12750_deeplane2_trimlowqual_trimadapterv2/gsnap.GM12750_deeplane2_trimlowqual_trimadapterv2.defv20120410.unique.bam,/gpfs/fs121/h/toung/aln/baseline/GM12750_deeplane3_trimlowqual_trimadapterv2/gsnap.GM12750_deeplane3_trimlowqual_trimadapterv2.defv20120410.unique.bam,/gpfs/fs121/h/toung/aln/baseline/GM12750_deeplane4_trimlowqual_trimadapterv2/gsnap.GM12750_deeplane4_trimlowqual_trimadapterv2.defv20120410.unique.bam,/gpfs/fs121/h/toung/aln/baseline/GM12750_deeplane5_trimlowqual_trimadapterv2/gsnap.GM12750_deeplane5_trimlowqual_trimadapterv2.defv20120410.unique.bam,/gpfs/fs121/h/toung/aln/baseline/GM12750_deeplane6_trimlowqual_trimadapterv2/gsnap.GM12750_deeplane6_trimlowqual_trimadapterv2.defv20120410.unique.bam,/gpfs/fs121/h/toung/aln/baseline/GM12750_deeplane7_trimlowqual_trimadapterv2/gsnap.GM12750_deeplane7_trimlowqual_trimadapterv2.defv20120410.unique.bam,/gpfs/fs121/h/toung/aln/baseline/GM12750_deeplane8_trimlowqual_trimadapterv2/gsnap.GM12750_deeplane8_trimlowqual_trimadapterv2.defv20120410.unique.bam,/gpfs/fs121/h/toung/aln/baseline/GM12750_onelane_trimlowqual_trimadapterv2/gsnap.GM12750_onelane_trimlowqual_trimadapterv2.defv20120410.unique.bam";
my $fasta_folder = "/gpfs/fs121/h/toung/work/2012-04-24_WholeFlowCell/step1-PoolRDD/step7-BlatReads/step1_ExtractRDDReads/test_fasta"; ### (optional)

my $options = GetOptions(
	"rddFile=s" => \$rddFile,
	"tag=s" => \$tag,
	"region=s" => \$region,
	"bam=s" => \$bam,
	"fasta_folder=s" => \$fasta_folder
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
my @bam = split(',',$bam);
foreach my $b (@bam) {
	(print STDERR "bam:\t$b\n") && -e $b or die "[STDERR]: bam doesn't exist $b\n";
}
if (defined $fasta_folder) {
	(print STDERR "fasta_folder:\t$fasta_folder\n") && -d $fasta_folder or die "[STDERR]: fasta folder doesn't exist\n";
}

my $fasta_file;
$fasta_file = $fasta_folder . "/" if defined $fasta_folder;
$fasta_file .= $tag . "_";
$fasta_file .= $regionPart . "_" . $regionStart . "_" . $regionEnd . "_" if defined $region;
print STDERR "fasta_file:\t$fasta_file\n";
my $fasta = Bio::SeqIO->new(-file => ">$fasta_file", -format => 'fasta');

open(FILE,$rddFile) or die "[STDERR]: can't open $rddFile: $!\n";

while(<FILE>) {
	chomp;

	next if (/^#/);
	
	my $lineNumber = $.;

	if (defined $region) {
		next if ($lineNumber < $regionStart);
		next if ($lineNumber > $regionEnd);
	}

	my @split = split('\t');
	my $chrom = $split[0];
	my $end = $split[2];
	
	my $region = $chrom . ":" . $end . "-" . $end;		
	foreach my $b (@bam) {
		my $command = "samtools view $b $region |";
		open(COMMAND,$command) or die "[STDERR]: can't open $command: $!\n";
	
		while(<COMMAND>) {
			chomp;
			
		}
		close(COMMAND);
	}
	

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
