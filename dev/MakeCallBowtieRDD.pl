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
my $exclude;
my $strand_specific;
my $adapter_only;
my $unique_aln_only;
my $unique_seq_only;
my $minqual;
my $min_rdd_level;
my $min_total_count;
my $max_size;
my $all_chrom;

my $options = GetOptions (
	"bam=s" => \$bam,
	"index=s" => \$index,
	"tag=s" => \$tag,
	"alnDB=s" => \$alnDB,
	"exclude=s" => \$exclude,
	"strand_specific=i" => \$strand_specific,
	"adapter_only=s" => \$adapter_only,
	"unique_aln_only=i" => \$unique_aln_only,
	"unique_seq_only=s" => \$unique_seq_only,
	"minqual=i" => \$minqual,
	"min_rdd_level=s" => \$min_rdd_level,
	"min_total_count=s" => \$min_total_count,
	"max_size=i" => \$max_size,
	"all_chrom=i" => \$all_chrom
);

### check that $BAM ends in '.bam' and is an absolute file name ################
my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
$bam_ext eq '.bam' or die "[STDERR]: '$bam' does not end in '.bam'\n";
substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$bam' is not absolute file name\n";
$bam = $bam_dir . $bam_name . ".unique" . $bam_ext if $unique_aln_only;
($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
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
	### get 'exclude' ready ################################################
	my @exclude_new;
	foreach my $EXCLUDE (split(',',$exclude)) {
		my $exclude_sub = $EXCLUDE;
		$exclude_sub =~ s/XXX/$CHROM/; ### substitute in the chromosome
		if(-e $exclude_sub) { push(@exclude_new,$exclude_sub); }
		else { print "\t\t[WARNING]: $exclude_sub doesn't exist\n"; }
	}
	my $exclude_new = join(',',@exclude_new);
	$exclude_new = 0 if $exclude_new eq "";
	########################################################################

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

		my $RUN_FILE = "CallBowtieRDD-$bam_name-$tag-$CHROM-$FILE->[0]-$FILE->[1].sh"; 
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

perl $HOME/dev/CallBowtieRDD.pl --bam $bam --index $index --region $REGION --tag $tag --alnDB $alnDB --exclude $exclude_new --strand_specific $strand_specific --adapter_only $adapter_only --unique_aln_only $unique_aln_only --unique_seq_only $unique_seq_only --minqual $minqual --min_rdd_level $min_rdd_level --min_total_count $min_total_count

END
	}
}

