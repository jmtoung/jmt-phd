#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Database;
use Bio::DB::Sam;
use File::Basename;
use Statistics::Descriptive;
use RegularList;
use PileupData;
use ComplementBase;

my $bam;
my $index;
my $chrom;
my $tag;
my $exclude;
my $min_rdd_level;
my $min_total_reads;
my $unique_only;
my $adapter_only;
my $strand_specific;

my $options = GetOptions("bam=s" => \$bam, "index=s" => \$index, "chrom=s" => \$chrom, "tag=s" => \$tag,
			"exclude=s" => \$exclude, "min_rdd_level=i" => \$min_rdd_level,
			"min_total_reads=i" => \$min_total_reads, "unique_only=i" => \$unique_only,
			"adapter_only=i" => \$adapter_only, "strand_specific=i" => \$strand_specific);

################################################################################
### This script outputs a list of candidate RDD sites from a bam file
################################################################################

### LOAD $bam AND MAKE SITES DIRECTORY #########################################
defined $index && -e $index or die "[STDERR]: missing $index: $!\n";
my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
$bam_ext eq '.bam' or die "[STDERR]: '$bam' does not end in '.bam'\n";
substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$bam' is not absolute file name\n";
### make a directory for the output that ends in ".sites"
my $output_dir = $bam_dir . $bam_name . ".$tag.sites";
unless(-d $output_dir) { mkdir $output_dir; }
chdir($output_dir);
### make output file
my $output_file = $output_dir . "/$bam_name.$tag.$chrom.sites";
open(OUTPUT,">$output_file") or die "[STDERR]: can't open $output_file: $!\n";
################################################################################

### LOAD $dna (all dna sites we want to exclude) ###############################
my @EXCLUDE;
if ($exclude) {
	foreach my $EXCLUDE (split(',',$exclude)) { push(@EXCLUDE,RegularList->new($EXCLUDE)); }
}
################################################################################

### OPEN RDD OUTPUT FILE #######################################################
my $pileup = "samtools view -h $bam '$chrom' | ";

if ($unique_only && !$adapter_only) {
	### awk '$1~/^@/ || $5==1'
	$pileup .= "awk '" . '$1~/^@/ || $5==1' . "' ";
} elsif ($adapter_only && !$unique_only) {
	### awk '$1~/^@/ || $1~/.*\|[0-9]+\|[0-9]+$/'
	$pileup .= "awk '" . '$1~/^@/ || $1~/.*\|[0-9]+\|[0-9]+$/' . "' ";
} elsif ($unique_only && $adapter_only) {
	### awk '$1~/^@/ || ($5==1 && $1~/.*\|[0-9]+\|[0-9]+$/)'
	$pileup .= "awk '" . '$1~/^@/ || ($5==1 && $1~/.*\|[0-9]+\|[0-9]+$/)' . "' ";
} else {
	### awk '$1~/^@/'
	$pileup .= "awk '" . '$1~/^@/ || $5>0' . "' ";
}

$pileup .= "| samtools view -bS - | samtools mpileup -f $index - |";
print "pileup command: $pileup\n";
open(PILEUP,$pileup) or die "[STDERR]: can't fork $pileup: $!\n";
################################################################################

################################################################################
while(<PILEUP>) {
	chomp;
	my @split = split('\t');

	my $CHROM = $split[0];
	my $POSITION = $split[1];
	
	foreach my $EXCLUDE (@EXCLUDE) {
		while($EXCLUDE->has_next && $EXCLUDE->get_column(3) < $POSITION) { $EXCLUDE->update_next_line; }
		next if ($EXCLUDE->get_column(2) <= $POSITION);
	}

	my $PILEUP = PileupData->new(\@split);

	my @STRANDS;
	if ($strand_specific) { push(@STRANDS,['+'],['-']); }
	else { push(@STRANDS,['+','-']); }

	foreach my $STRANDS (@STRANDS) {

		my $REF_BASE = $split[2];

		my $TOTAL_COUNT = $PILEUP->get_base_count($STRANDS,['A','C','G','T']);
		goto START if $TOTAL_COUNT == 0;
		if ($min_total_reads) { goto START if $TOTAL_COUNT < $min_total_reads; }

		foreach my $BASE ('A','C','T','G') {

			if ($STRANDS->[0] eq '-') { goto START if $BASE eq complement_base($REF_BASE); }
			else { goto START if $BASE eq $REF_BASE; }

			my $BASE_COUNT = $PILEUP->get_base_count($STRANDS,[$BASE]);
			goto START if $BASE_COUNT == 0;

			my $RDD_LEVEL = $BASE_COUNT / $TOTAL_COUNT;

			if($min_rdd_level) { goto START if $RDD_LEVEL < $min_rdd_level; }

			print OUTPUT $CHROM, "\t", $POSITION, "\n";
			goto START;
		}
	}
	START:
}
