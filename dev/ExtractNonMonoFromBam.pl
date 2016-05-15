#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Database;
use File::Basename;
use PileupData;
use ComplementBase;

my $bam;
my $index;
my $chrom;
my $minqual;
my $tag;

my $options = GetOptions(
	"bam=s" => \$bam,
	"index=s" => \$index,
	"chrom=s" => \$chrom,
	"minqual=i" => \$minqual,
	"tag=s" => \$tag
);

################################################################################
### This script outputs a list of candidate non-monomorphic sites from a bam file
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

### OPEN RDD OUTPUT FILE #######################################################
print "run options: perl ~/dev/ExtractNonMonoFromBam.pl --bam $bam --index $index --chrom $chrom --tag $tag --minqual $minqual\n\n";
my $pileup = "samtools mpileup -r $chrom -f $index $bam |";
print "pileup command: $pileup\n\n";
open(PILEUP,$pileup) or die "[STDERR]: can't fork $pileup: $!\n";
################################################################################

################################################################################
while(<PILEUP>) {
	chomp;
	my @split = split('\t');
	my $CHROM = $split[0];
	my $POSITION = $split[1];
	my $REF_BASE = $split[2];

	my $PILEUP = PileupData->new(\@split,$minqual);
	my $REF_BASE_COUNT = $PILEUP->get_base_count(['+','-'],[$REF_BASE]);
	my @NONREF_BASE;
	foreach my $BASE ('A','C','G','T') { push(@NONREF_BASE,$BASE) if $BASE ne $REF_BASE; }
	my $NONREF_BASE_COUNT = $PILEUP->get_base_count(['+','-'],\@NONREF_BASE);
	my $DELETION_COUNT = $PILEUP->get_base_count(['+','-'],['X']);
	my $INSERTION_COUNT = $PILEUP->get_ins_count(['+','-']);
	print OUTPUT $CHROM, "\t", $POSITION, "\t", $REF_BASE_COUNT, "\t", $NONREF_BASE_COUNT, "\t", $DELETION_COUNT, "\t", $INSERTION_COUNT, "\n";
}

my @split = split('\/',$bam);
my @entry = ($split[-1],$chrom);
my $primary_key = $split[-1] . "-" . $chrom;
Database->new("/ifs/h/toung/rdd/1000_genomes/pilot_data/finished_pileup_minqual20_db")->add2db(\@entry,0,$primary_key);
