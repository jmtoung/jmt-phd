#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib '/gpfs/fs121/h/toung/dev', '/gpfs/fs121/h/toung/lib/perl5/share/perl5';
use Bio::DB::Sam;
use File::Basename;
use Statistics::Descriptive;
use CalculateRDDStats_v20121203;
use Data::Dumper;

umask 0007;

my $sites; ### REQUIRED
my $region; ### optional
my $bam; ### full path please
my $index; ### REQUIRED
my $strand_specific; ### REQUIRED
my $unique_aln_only; ### REQUIRED
my $unique_seq_only; ### REQUIRED
my $adapter_only; ### REQUIRED
my $minqual = 20; ### REQUIRED
my $alnID; ### required
my $combine_trim = 0;
my $edit_dist_filter = 0;

$|++;
my $result = GetOptions(
	"sites=s" => \$sites,
	"region=s" => \$region,
	"bam=s" => \$bam,
	"index=s" => \$index,
	"strand_specific=i" => \$strand_specific,
	"unique_aln_only=i" => \$unique_aln_only,
	"unique_seq_only=s" => \$unique_seq_only,
	"adapter_only=s" => \$adapter_only, 
	"minqual=i" => \$minqual,
	"alnID=s" => \$alnID,
	"combine_trim=i" => \$combine_trim,
	"edit_dist_filter=i" => \$edit_dist_filter,
);

################################################################################
### This script extracts RDD info for a bam file given a list of sites & RDD types
################################################################################

### PRINT OPTIONS ##############################################################
print STDERR "sites:\t$sites\n";
print STDERR "region:\t$region\n" if defined $region;
(-e $index && print STDERR "index:\t$index\n") or die "[STDERR]: index $index not defined\n";
($strand_specific == 0 || $strand_specific == 1) && print STDERR "strand_specific:\t$strand_specific\n" or die "[STDERR]: not defined strand_specific\n";
($unique_aln_only == 0 || $unique_aln_only == 1) && print STDERR "unique_aln_only:\t$unique_aln_only\n" or die "[STDERR]: not defined unique_aln_only\n";
($unique_seq_only eq "0,1" || $unique_seq_only eq "0" || $unique_seq_only eq "1") && print STDERR "unique_seq_only:\t$unique_seq_only\n" or die "[STDERR]: not defined unique_seq_only\n";
($adapter_only eq "0,1" || $adapter_only eq "0" || $adapter_only eq "1") && print STDERR "adapter_only:\t$adapter_only\n" or die "[STDERR]: not defined adapter_only\n";
defined $minqual && $minqual =~ /^[0-9]+$/ && print STDERR "minqual:\t$minqual\n" or die "[STDERR]: not defined minqual\n";
($combine_trim == 0 || $combine_trim == 1) && print STDERR "combine_trim:\t$combine_trim\n" or die "[STDERR]: combine trim not defined\n";
($edit_dist_filter == 0 || $edit_dist_filter == 1) && print STDERR "edit_dist_filter:\t$edit_dist_filter\n" or die "[STDERR]: edit_dist_filter not defined\n"; ### added 7/13/2012 2pm
defined $alnID && print STDERR "alnID:\t$alnID\n" or die "[STDERR]: alnID not defined\n";
################################################################################

### LOAD BAM & INDEX ###########################################################
my @bam = split(',',$bam);
my @BAM; ### bam objects
foreach my $b (@bam) {
	print STDERR "bam:\t$b\n";
	if ($unique_aln_only) {
		($b =~ /unique/ || $b =~ /primary/) or die "[STDERR]: specified unique_aln_only but not primary or unique\n";
	}
	my $BAM = Bio::DB::Sam->new(-bam => $b,-fasta => $index,-autoindex => 1);
	push(@BAM,$BAM);
}
print STDERR "numBamFiles:\t", scalar(@bam), "\n";
################################################################################

### CHECK IF REGION IS DEFINED #################################################
my ($region_index,$region_start,$region_end);
if (defined $region) {
	if ($region =~ /([0-9]+):([0-9]+)-([0-9]+)/) { $region_index = $1; $region_start = $2; $region_end = $3; } 
	else { die "improper region $region\n"; }
}
################################################################################

### LOAD NONREF/RDD SITES ######################################################
open(SITES,$sites) or die "[STDERR]: can't open $sites: $!\n";
while(<SITES>) {
	chomp;
	
	if (defined $region) { 
		last if $. > $region_end;
		next unless $. >= $region_start;
	}

	my @split = split('\t');

	next if ($_ =~ /^#/ || $_ =~ /^chrom/);

	my $CHROM = $split[0];
	my $START = $split[1];
	my $END = $split[2];
	my $STRAND = $split[3];
	my $REF_BASE = $split[4];
	my $RDD_BASE = $split[5];

	my @strands;
	if ($strand_specific) {
		$STRAND eq '+' or $STRAND eq '-' or die "[STDERR]: strand specific but strand is $STRAND\n";
		push(@strands,[$STRAND]);
	} else {
		push(@strands,['+','-']);
		if ($STRAND eq '-') {
			$REF_BASE =~ tr/ACGT/TGCA/;
			$RDD_BASE =~ tr/ACGT/TGCA/;
		}
	}

	my $RDD_OBJ = makeRddObject(\@BAM,$CHROM,$END,$minqual,$combine_trim,$edit_dist_filter);	

	foreach my $STRANDS (@strands) {

		foreach my $ADAPTER_ONLY (split(',',$adapter_only)) {

			foreach my $UNIQUE_SEQ_ONLY (split(',',$unique_seq_only)) {

				my @RESULTS;				
				push(@RESULTS,$CHROM,$START,$END,join(",",@{$STRANDS}),$REF_BASE,$RDD_BASE);
				push(@RESULTS,$unique_aln_only,$UNIQUE_SEQ_ONLY,$ADAPTER_ONLY);

				### GET ALL BASES
				my $TOTAL_BASES = getBases($RDD_OBJ,$STRANDS);

				### DO THE FOLLOWING FOR TWO GROUPS OF BASES
				foreach my $GROUP ($TOTAL_BASES,[$RDD_BASE]) {

					my $STATS = getRddStats($RDD_OBJ,$GROUP,$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY);
					foreach my $STAT (@{$STATS}) {
						if (ref($STAT) eq '') {
							push(@RESULTS,$STAT);
						} elsif (ref($STAT) eq 'ARRAY') {
##							my $CALC = calcStats($STAT);
##							push(@RESULTS,join(";",@{$CALC}));
							push(@RESULTS,join(";",@{calcStats($STAT)}));
						} else {
							die "[STDERR]: weird referenced object!\n";
						}
					}
				}
				
				push(@RESULTS,$alnID);
				
				print join("\t",@RESULTS), "\n";
			}
		}
	}
}
close(SITES);


################################################################################

sub calcStats {
	my $ARRAY = shift;

	if (scalar(@{$ARRAY}) != 0) {
		my $STAT = Statistics::Descriptive::Full->new();
		$STAT->add_data($ARRAY);
		my $stdev = 'NA';
		if (scalar(@{$ARRAY}) > 1) {
			$stdev = $STAT->standard_deviation;
		}
		return [$STAT->mean(), $stdev, $STAT->max(), $STAT->min()];
	} else {
		return ['NA','NA','NA','NA'];
	}
}
