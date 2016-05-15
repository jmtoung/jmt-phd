#!/usr/bin/perl -w

use strict;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Getopt::Long;
use File::Basename;
use PileupData;
use ComplementBase;

my $list;
my $index;
my $bam;
my $tag;

my $options = GetOptions("list=s" => \$list, "index=s" => \$index, "bam=s" => \$bam, "tag=s" => \$tag);

################################################################################
### This script takes a list of candidate RDD sites, and a list
### of bam files (DNA) and calculates pileup on that list, filtering out sites 
### where the RDD allele is present in the DNA sequencing information
################################################################################

### check if index exists & file exists
-e $index or die "[STDERR]: $index doesn't exist: $!\n";
-e $list or die "[STDERR]: $list doesn't exist: $!\n";
my ($list_name, $list_dir, $list_ext) = fileparse($list,'\.sites');
$list_ext eq '.sites' or die "[STDERR]: '$list' does not end in '.sites'\n";
### go to directory of list file & make a folder with the $tag name
chdir($list_dir) or die "[STDERR]: can't go to $list_dir: $!\n";
mkdir $tag unless(-d $tag);
chdir($tag) or die "[STDERR]: can't go to $tag: $!\n";
open(OUTPUT,">$list_name") or die "[STDERR]: can't open $list_name: $!\n";
select(OUTPUT);
################################################################################

### read in the bam file and check that files exist 
open(BAM,$bam) or die "[STDERR]: can't open file: $!\n";
chomp(my $BAM = <BAM>); ### read in first line
my @BAM = split(' ',$BAM);
my @BAM_new;
foreach my $file (@BAM) {
	push(@BAM_new,$file) if -e $file;
}
@BAM = @BAM_new;
my $BAM_list = join(" ",@BAM);
################################################################################

my $LIST = RegularList->new($list);

### get chromosomes
my $pileup = "samtools mpileup -f $index -l $list $BAM_list |";
open(PILEUP,$pileup) or die "[STDERR]: can't fork $pileup: $!\n";
while(<PILEUP>) {
	chomp;
	my @split = split('\t');
	
	my $CHROM = $split[0];
	my $POSITION = $split[1];
	my $REF_BASE = $split[2];

	while($LIST->has_next && $LIST->get_column(1) < $POSITION) {
		print $LIST->get_column(0), "\t", $LIST->get_column(1), "\n";
		$LIST->update_next_line;
	} 

	if ($LIST->has_next && $LIST->get_column(1) == $POSITION) {
		$DB::single = 1;
		my $PILEUP = PileupData->new(\@split);
		my @NONREF_BASE;
		foreach my $BASE ('A','C','G','T') { push(@NONREF_BASE,$BASE) if $BASE ne $REF_BASE; }
		my $TOTAL_COUNT = $PILEUP->get_base_count(['+','-'],\@NONREF_BASE);
		
	} 

}

while($LIST->has_next) {
	print $LIST->get_column(0), "\t", $LIST->get_column(1), "\n";
}
################################################################################
