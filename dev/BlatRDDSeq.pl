#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use File::Basename;
use PullGenomeSequence;
use Bio::SeqIO;
use Bio::Seq;

my $HOME;
my $files;
my $index;
my $tag;
my $flanking_distance;

my $options = GetOptions(
	"HOME=s" => \$HOME,
	"files=s" => \$files,
	"index=s" => \$index,
	"tag=s" => \$tag,
	"flanking_distance=s" => \$flanking_distance
);

($HOME eq '/ifs/h/toung' or $HOME eq '/home/jmtoung/Lab') or die "[STDERR]: defined HOME\n";
print "HOME:\t$HOME\n";

my @files = split(',',$files);
$files = join(" ",@files);
my @flanking_distance = split(',',$flanking_distance);

foreach my $file (@files) { 
	-e $file or die "[STDERR]: $file doesn't exist\n";	
	print "file:\t$file\n"; 
}

print "index:\t$index\n";
print "tag:\t$tag\n";

foreach my $fl_dist (@flanking_distance) { $fl_dist =~ /[0-9]+/ or die "[STDERR]: flanking distance invalid\n"; }
print "flanking_distance:\t$flanking_distance\n";

### make fasta files with sequences for each flanking distance
my %fasta;
foreach my $fl_dist (split(',',$flanking_distance)) {
	my $fasta_file = "BlatRDDSeq-$tag-$fl_dist.fasta";
	print "fasta_file:\t$fasta_file\n";
	$fasta{$fl_dist}{'fasta_file'} = $fasta_file;
	$fasta{$fl_dist}{'fasta_obj'} = Bio::SeqIO->new(-file => ">$fasta{$fl_dist}{'fasta_file'}", -format => 'fasta');
}

### read the unique sites from the files and create DNA and RNA sequences
my $command = "cat $files | sort +0 -1 +1n -2 +2 -3 +3 -4 +5 -6 |";
print "cat_command:\t$command\n";

open(COMMAND,$command) or die "[STDERR]: can't fork $command: $!\n";
my %current; ### hold current RDD
while(<COMMAND>) {
	chomp;
	my @split = split('\t');
	$split[2] eq '+,-' or die "[STDERR]: strand doesn't equal '+,-'\n";
	
	if (is_equal(\@split,\%current)) { 
		add_to_current(\@split,\%current);
	} else {
		make_sequences(\%current,$flanking_distance,\%fasta,$index);
		%current = ();
		add_to_current(\@split,\%current);
	}
}
make_sequences(\%current,$flanking_distance,\%fasta,$index);

## RUN BLAT ON EACH OF THE FASTA FILES ########################################
my %blat_results;
foreach my $fl_dist (split(',',$flanking_distance)) {
	my $blat_output = $fasta{$fl_dist}{'blat'} = "BlatRDDSeq-$tag-$fl_dist.psl";

	my $blat = "blat";
	if ($HOME eq '/ifs/h/toung') { $blat = $HOME . "/bin/" . $blat; }

	print "\nrunning blat ($blat)...\n";
	print "blat_command:\ttime $blat -stepSize=5 -minIdentity=95 -repMatch=2253 -noHead $index $fasta{$fl_dist}{'fasta_file'} $blat_output\n";
	
	system("time $blat -stepSize=5 -minIdentity=95 -repMatch=2253 -noHead $index $fasta{$fl_dist}{'fasta_file'} $blat_output");

	### parse results from blat output
	open(BLAT,$blat_output) or die "[STDERR]: can't open $blat_output: $!\n";

	while(<BLAT>) {
		chomp;
		my @split = split('\t');
		
		### require total number of bases that match equal the size of the query sequence
		next unless ($split[0] + $split[2] >= ($split[10] - 2));
		
		### require number of inserts in query and number of bases inserted in query both be 0
		next unless ($split[4] == 0 && $split[5] == 0);
		
		### require number of inserts in target and number of bases inserted in target both be 0
		next unless ($split[6] == 0 && $split[7] == 0);
		
		### require block count be 1 (no gaps)
		next unless $split[17] == 1;

		### query info
		my ($chrom,$position,$flank_dist,$ref_base,$rdd_base,$mol) = split(';',$split[9]);
		my $query_start = $position - $flank_dist - 1;
		my $query_end = $position + $flank_dist;
		
		### target info
		my ($target_chrom,$target_start,$target_end) = ($split[13],$split[15],$split[16]);
	
		### next if alignment is in the same location as the RDD.
		next if ($chrom eq $target_chrom && $query_start == $target_start && $query_end == $target_end);
		
		$blat_results{$chrom}{$position}{$mol}{$flank_dist}++;
	}
	close(BLAT);
}


### open each file
foreach my $file (@files) {

	open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";
	my ($file_name,$file_dir,$file_ext) = fileparse($file,'\.txt');
	my $output = $file_dir . $file_name . "_blat" . $file_ext;
	open(OUTPUT,">".$output) or die "[STDERR]: can't open $output: $!\n";
	
	while(<FILE>) {
		chomp;
		my @split = split('\t');
		my $chrom = $split[0];
		my $position = $split[1];
		my $ref_base = $split[3];
		my $rdd_base = $split[5];
		
		foreach my $mol ('DNA','RNA') {
			foreach my $fl_dist (@flanking_distance) {
				if (exists $blat_results{$chrom}{$position}{$mol}{$fl_dist}) {
					push(@split,$blat_results{$chrom}{$position}{$mol}{$fl_dist});
				} else {
					push(@split,0);
				}
			}
		}
		
		print OUTPUT join("\t",@split), "\n";
	}
	close(FILE);
	close(OUTPUT);
}

sub is_equal {
	my ($split,$current) = @_;
	
	### if hash is empty (like for first one), return equal
	return 1 if !exists $current->{'chrom'};
	
	return 0 unless $current->{'chrom'} eq $split->[0];
	return 0 unless $current->{'position'} eq $split->[1];
	return 0 unless $current->{'strands'} eq $split->[2];
	return 0 unless $current->{'ref_base'} eq $split->[3];
	
}

sub add_to_current {
	my ($split,$current) = @_;
	
	$current->{'chrom'} = $split->[0];
	$current->{'position'} = $split->[1];
	$current->{'strands'} = $split->[2];
	$current->{'ref_base'} = $split->[3];
	$current->{'rdd_base'}{$split->[5]}++;
}

sub make_sequences {
	my ($current,$flanking_distance,$fasta,$index) = @_;
	
	foreach my $fl_dist (split(',',$flanking_distance)) {
		my $display_id = join(";",$current->{'chrom'},$current->{'position'},$fl_dist,$current->{'ref_base'});
		my $region = $current->{'chrom'} . ":" . ($current->{'position'} - $fl_dist) . "-" . ($current->{'position'} + $fl_dist);
		my $sequence = pull_sequence($region,$index);
		### sequence may not be defined if at ends of chromosome
		next unless $sequence;

		### DNA SEQUENCE			
		my $dna_obj = Bio::Seq->new(-seq => $sequence, -display_id => $display_id . ";" . $current->{'ref_base'} . ";DNA", -desc => "", -alphabet => "dna");
		$fasta->{$fl_dist}{'fasta_obj'}->write_seq($dna_obj);

		### RNA SEQUENCE
		my @sequence = split('',$sequence);
		foreach my $rdd_base (sort keys %{$current->{'rdd_base'}}) {
			my @RNA_sequence = @sequence;
			$RNA_sequence[$fl_dist] = $rdd_base;
			my $RNA_sequence = join('',@RNA_sequence);
		
			my $rna_obj = Bio::Seq->new(-seq => $RNA_sequence, -display_id => $display_id . ";$rdd_base;RNA", -desc => "", -alphabet => "dna");
			$fasta->{$fl_dist}{'fasta_obj'}->write_seq($rna_obj);
		}
	}
}
