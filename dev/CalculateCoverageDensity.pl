#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Database;
use File::Basename;
use PileupData;

my $bam;
my $index;
my $region;
my $tag;
my $alnDB;
my $unique_aln_only;
my $minqual;

my $result = GetOptions(
	"bam=s" => \$bam,
	"index=s" => \$index,
	"region=s" => \$region, 
	"tag=s" => \$tag,
	"alnDB=s" => \$alnDB,
	"unique_aln_only=i" => \$unique_aln_only,
	"minqual=i" => \$minqual,
);

################################################################################
### This script calculates the coverage density for a bam file
################################################################################

### LOOK UP bamID ##############################################################
my $bamID = Database->new($alnDB)->lookup(1,$bam,0);
defined $bamID or die "[STDERR]: bamID for '$bam' not defined in '$alnDB'\n";
################################################################################

### LOAD BAM ###################################################################
my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
$bam_ext eq '.bam' or die "[STDERR]: '$bam' does not end in '.bam'\n";
substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$bam' is not absolute file name\n";
$bam_name .= ".unique" if $unique_aln_only;
$bam = $bam_dir . $bam_name . $bam_ext;
defined $bam && -e $bam or die "[STDERR]: missing $bam: $!\n";
defined $index && -e $index or die "[STDERR]: missing $index: $!\n";
################################################################################

### MAKE OUTPUT DIRECTORY & FILES ##############################################
my $cov_dir = $bam_dir . $bam_name;
$cov_dir =~ s/\/aln\//\/cov\//;
unless(-d $cov_dir) { system("mkdir -p $cov_dir"); }
my $region_sub = $region;
$region_sub =~ s/:/\./;
my $output = $cov_dir . "/$bam_name.$tag.$region_sub.cov";
open(OUTPUT,">".$output) or die "[STDERR]: can't open $output: $!\n";
select(OUTPUT);
################################################################################

################################################################################
my $pileup;
### if $region is a file, then we want -l option, else want -r option
if (-e $region) { $pileup = "samtools mpileup -f $index -l $region $bam |"; }
else { $pileup = "samtools mpileup -f $index -r $region $bam |"; }
print STDOUT "pileup command: $pileup\n";
open(PILEUP,$pileup) or die "[STDERR]: can't open $pileup: $!\n";
my %COV;
while(<PILEUP>) {
	chomp;
	my @split = split('\t');

	my $CHROM = $split[0];
	my $POSITION = $split[1];
	
	my $RNA_pileup = PileupData->new(\@split,$minqual);
	next unless $RNA_pileup->is_def;

	my @STRANDS = (['+','-']);

	foreach my $STRANDS (@STRANDS) {

		my $TOTAL_COUNT = $RNA_pileup->get_base_count($STRANDS,['A','C','G','T']);
		next if $TOTAL_COUNT == 0;
		
		$COV{$TOTAL_COUNT}++;
	}
}

foreach my $TOTAL_COUNT (sort {$a <=> $b} keys %COV) {
	print $TOTAL_COUNT, "\t", $COV{$TOTAL_COUNT}, "\t", $bamID, "\n";
}
