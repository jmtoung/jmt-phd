#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use Cwd;

my $bam;
my $index;
my $disableBAQ = 1;

my $result = GetOptions(
	"bam=s" => \$bam,
	"index=s" => \$index,
	"disableBAQ=i" => \$disableBAQ
);

print STDERR "bam:\t$bam\n";
print STDERR "index:\t$index\n";
print STDERR "disableBAQ:\t$disableBAQ\n";

################################################################################
### This script takes a bam file and outputs pileup on all of its chromosomes
################################################################################

### check if bam exists ########################################################
my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.[a-zA-Z0-9]+$');
$bam_ext eq '.bam' or die "[STDERR]: '$bam' does not end in '.bam'\n";
substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$bam' is not absolute file name\n";
-e $bam or die "[STDERR]: $bam doesn't exist\n";
################################################################################

### make pileup directory ######################################################
my $pileup_dir = $bam_dir . $bam_name . ".pileup";
unless (-e $pileup_dir) { mkdir $pileup_dir; }
################################################################################

### check index exists #########################################################
-e $index or die "[STDERR]: '$index' does not exist\n";
################################################################################

### pileup #####################################################################
my $pileup = "samtools mpileup ";
$pileup .= "-B " if $disableBAQ;
$pileup .= "-f $index $bam |";
print STDERR $pileup, "\n";
open(PILEUP,$pileup) or die "[STDERR]: can't fork $pileup: $!\n";

my %output;
while(<PILEUP>) {
	chomp;
	my @split = split('\t');
	my $chrom = $split[0];

	unless (exists $output{$chrom}) {
		my $output = $pileup_dir . "/" . $bam_name . ".$chrom.pileup";
		open($output{$chrom},">".$output) or die "[STDERR]: can't open $output: $!\n";
	}
	print {$output{$chrom}} $_, "\n";
}
