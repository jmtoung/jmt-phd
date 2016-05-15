#!/usr/bin/perl

use strict;
use lib '/ifs/h/toung/dev', '/home/jmtoung/Lab/dev';
use List::Util qw(max min);
use Database;
use File::Basename;
use Getopt::Long;

my $HOME = "/ifs/h/toung";

my $bam;
my $index;
my $tag;
my $alnDB;
my $unique_aln_only;
my $minqual;
my $max_size;
my $all_chrom;

my $options = GetOptions (
	"bam=s" => \$bam,
	"index=s" => \$index,
	"tag=s" => \$tag,
	"alnDB=s" => \$alnDB,
	"unique_aln_only=i" => \$unique_aln_only,
	"minqual=i" => \$minqual,
	"max_size=i" => \$max_size,
	"all_chrom=i" => \$all_chrom
);

### check that $BAM ends in '.bam' and is an absolute file name ################
my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
$bam_ext eq '.bam' or die "[STDERR]: '$bam' does not end in '.bam'\n";
substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$bam' is not absolute file name\n";
-e $bam or die "[STDERR]: $bam doesn't exist: $!\n";
################################################################################

### get chromosomes in $bam ####################################################
my %CHROM_BAM;
my $bam_bai = "samtools idxstats $bam |";
open(BAM_BAI,$bam_bai) or die "[STDERR]: can't fork '$bam_bai': $!\n";
while(<BAM_BAI>) {
	chomp;
	my @split = split('\t');
	$CHROM_BAM{$split[0]} = $split[1]; ### $CHROM_BAM{'chrom1'} = size_of_chrom1
}
################################################################################

### get chromosomes in $index ##################################################
my %CHROM_INDEX;
my $index_fai = $index . ".fai";
open(INDEX_FAI,$index_fai) or die "[STDERR]: can't open fai_index: $index_fai: $!\n";
while(<INDEX_FAI>) {
	chomp;
	my @split = split('\t');
	$CHROM_INDEX{$split[0]} = $split[1];
}
################################################################################

print ">>> Looping through chromosomes <<<\n";
foreach my $CHROM (keys %CHROM_BAM) {

	unless($all_chrom) { next unless $CHROM =~ /^(chr)?([0-9XY]+$|^MT$)/; }
	next unless exists $CHROM_INDEX{$CHROM}; ### make sure chrom in index

	print "\t...working on chromosome '$CHROM'...\n";

	### divy up each chromosome according to max_size & length of chrom ####
	my $START = 1;
	my $END = min($CHROM_BAM{$CHROM}, $max_size);
	my @FILES;
	do {
		my @REGION = ($START,$END);
		push(@FILES,\@REGION);
		$START = ++$END;
		$END = min($CHROM_BAM{$CHROM}, $START + $max_size - 1);
	} until ($START > $CHROM_BAM{$CHROM});
	########################################################################

	### print out run scripts for each file ################################
	foreach my $FILE (@FILES) {

		my $RUN_FILE = "CalculateCoverageDensity-$bam_name-$tag-$CHROM-$FILE->[0]-$FILE->[1].sh"; 
		my $REGION = $CHROM . ":" . $FILE->[0] . '-' . $FILE->[1];
		open(RUN_FILE,">$RUN_FILE") or die "[STDERR]: can't open $RUN_FILE: $!\n";

print RUN_FILE <<"END";
#\$ -cwd ### use current directory
#\$ -S /bin/bash ### program to execute script
#\$ -M toung\@mail.med.upenn.edu ### email address
##\$ -m ea ### mail is to be sent at abort and end time
#\$ -j y ### combine stdout and stderr
#\$ -t 1-1 ### array job #
#\$ -V ### use current environment variables
##\$ -pe DJ 12 ### parallel threads
#\$ -l mem_free=4G ### request memory

perl $HOME/dev/CalculateCoverageDensity.pl --bam $bam --index $index --region $REGION --tag $tag --alnDB $alnDB --unique_aln_only $unique_aln_only --minqual $minqual

END
	}
}

