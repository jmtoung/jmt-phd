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
use CalculateRDDStats;
use ComplementBase;

my $bam;
my $index;
my $region;
my $tag;
my $alnDB;
my $exclude;
my $strand_specific; ### 1 (yes) or 0 (no)
my $adapter_only; ### 1 means only adapter reads; 0 means everything; 0,1 means do both
my $unique_aln_only;  ### 1 means only unique sequences; 0 means everything
my $unique_seq_only; ### 1 means collapse unique seq to one; 0 means don't; 0,1 means do both
my $minqual;
my $min_rdd_level;
my $min_total_count;

my $result = GetOptions(
	"bam=s" => \$bam,
	"index=s" => \$index,
	"region=s" => \$region, 
	"tag=s" => \$tag,
	"alnDB=s" => \$alnDB,
	"exclude=s" => \$exclude,
	"strand_specific=i" => \$strand_specific,
	"adapter_only=s" => \$adapter_only, 
	"unique_aln_only=i" => \$unique_aln_only,
	"unique_seq_only=s" => \$unique_seq_only,
	"minqual=i" => \$minqual,
	"min_rdd_level=s" => \$min_rdd_level,
	"min_total_count=s" => \$min_total_count
);

################################################################################
### This script calls RDDs for Bowtie alignments
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
my $BAM = Bio::DB::Sam->new(-bam => $bam,-fasta => $index); ### will exit if bam doesn't exist.
defined $index && -e $index or die "[STDERR]: missing $index: $!\n";
################################################################################

### MAKE OUTPUT DIRECTORY & FILES ##############################################
my $rdd_dir = $bam_dir . $bam_name;
$rdd_dir =~ s/\/aln\//\/rdd\//;
unless(-d $rdd_dir) { system("mkdir -p $rdd_dir"); }
my $region_sub = $region;
$region_sub =~ s/:/\./;
my $output = $rdd_dir . "/$bam_name.$tag.$region_sub.rdd";
open(OUTPUT,">".$output) or die "[STDERR]: can't open $output: $!\n";
select(OUTPUT);
$|++;
################################################################################

### LOAD sites we want to exclude) #############################################
my @EXCLUDE;
if ($exclude) { 
	foreach my $EXCLUDE (split(',',$exclude)) { push(@EXCLUDE,RegularList->new($EXCLUDE)); }
}
################################################################################

my ($region_chrom,$region_bound) = split(':',$region);
my ($region_start,$region_end) = split('-',$region_bound);

################################################################################
my $pileup;
### if $region is a file, then we want -l option, else want -r option
#if (-e $region) { $pileup = "samtools mpileup -f $index -l $region $bam |"; }
#else { $pileup = "samtools mpileup -f $index -r $region $bam |"; }
$pileup = "samtools mpileup -f $index $bam |"; 
print STDOUT "pileup command: $pileup\n";
open(PILEUP,$pileup) or die "[STDERR]: can't open $pileup: $!\n";
while(<PILEUP>) {
	chomp;
	my @split = split('\t');

	my $CHROM = $split[0];
	my $POSITION = $split[1];

	next unless ($CHROM eq $region_chrom && $POSITION >= $region_start && $POSITION <= $region_end);

	foreach my $EXCLUDE (@EXCLUDE) {
		while($EXCLUDE->has_next && $EXCLUDE->get_column(3) < $POSITION) { $EXCLUDE->update_next_line; }
		goto NEXT if ($EXCLUDE->has_next && $EXCLUDE->get_column(2) <= $POSITION);
	}
	
	my $RNA_pileup = PileupData->new(\@split,$minqual);
	goto NEXT unless $RNA_pileup->is_def;

	my @STRANDS;
	if ($strand_specific) { push(@STRANDS,['+'],['-']); }
	else { push(@STRANDS,['+','-']); }

	my $REF_BASE = $split[2];
	my $STATS_OBJ;

	foreach my $STRANDS (@STRANDS) {

		my $TOTAL_COUNT = $RNA_pileup->get_base_count($STRANDS,['A','C','G','T']);
		next if $TOTAL_COUNT == 0;
		if ($min_total_count) { next if $TOTAL_COUNT < $min_total_count; }

		foreach my $BASE ('A','C','T','G') {
			my @RESULTS;

			if ($STRANDS->[0] eq '-') { next if $BASE eq complement_base($REF_BASE); }
			else { next if $BASE eq $REF_BASE; }

			my $BASE_COUNT = $RNA_pileup->get_base_count($STRANDS,[$BASE]);
			next if $BASE_COUNT == 0;
			if ($min_rdd_level) { next if $BASE_COUNT/$TOTAL_COUNT < $min_rdd_level; }

			### CHROM & POSITION
			push(@RESULTS,$CHROM,$POSITION);

			### STRAND
			if (scalar(@{$STRANDS}) == 2) { push(@RESULTS,'*'); }
			elsif ($STRANDS->[0] eq '+') { push(@RESULTS,'+'); }
			else { push(@RESULTS,'-'); }
			
			### REF & RDD BASE
			if ($STRANDS->[0] eq '-') { push(@RESULTS,complement_base($REF_BASE),$BASE); }
			else { push(@RESULTS,$REF_BASE,$BASE); }

			$STATS_OBJ = make_RDD_stats_object($BAM,$CHROM,$POSITION,$REF_BASE,$minqual) if scalar(keys %{$STATS_OBJ}) == 0;

			foreach my $ADAPTER_ONLY (split(',',$adapter_only)) {
				foreach my $UNIQUE_SEQ_ONLY (split(',',$unique_seq_only)) {			

					### STATS FOR ALL BASES
					my ($COUNT,$QUAL,$MAPQUAL,$POS5,$POS3) = get_RDD_stats($STATS_OBJ,['A','C','G','T'],$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY);

					### STATS FOR RDD BASE
					my ($RDD_COUNT,$RDD_QUAL,$RDD_MAPQUAL,$RDD_POS5,$RDD_POS3) = get_RDD_stats($STATS_OBJ,[$BASE],$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY);

					next if $COUNT == 0 && $RDD_COUNT == 0;

					### PRINT RESULTS
					print join("\t",@RESULTS), "\t", $unique_aln_only, "\t", $UNIQUE_SEQ_ONLY, "\t", $ADAPTER_ONLY, "\t";

					print $COUNT, "\t";
					print "{";
					print join(",",@{calc_stats($QUAL)}) if scalar(@{$QUAL}) != 0;
					print "}", "\t";
					print "{";
					print join(",",@{calc_stats($MAPQUAL)}) if scalar(@{$MAPQUAL}) != 0;
					print "}", "\t";
					print "{";
					print join(",",@{calc_stats($POS5)}) if scalar(@{$POS5}) != 0;
					print "}", "\t";
					print "{";
					print join(",",@{calc_stats($POS3)}) if scalar(@{$POS3}) != 0;
					print "}", "\t";

					print $RDD_COUNT, "\t";
					print "{";
					print join(",",@{calc_stats($RDD_QUAL)}) if scalar(@{$RDD_QUAL}) != 0;
					print "}", "\t";
					print "{";
					print join(",",@{calc_stats($RDD_MAPQUAL)}) if scalar(@{$RDD_MAPQUAL}) != 0;
					print "}", "\t";
					print "{";
					print join(",",@{calc_stats($RDD_POS5)}) if scalar(@{$RDD_POS5}) != 0;
					print "}", "\t";
					print "{";
					print join(",",@{calc_stats($RDD_POS3)}) if scalar(@{$RDD_POS3}) != 0;
					print "}", "\t";
					
					print $bamID, "\n";
				}
			}
		}
	}
	NEXT:
}

sub calc_stats {
	my $ARRAY = shift;

	my $STAT = Statistics::Descriptive::Full->new();
	$STAT->add_data($ARRAY);
	return [$STAT->mean(), $STAT->standard_deviation, $STAT->max(), $STAT->min()];
}
