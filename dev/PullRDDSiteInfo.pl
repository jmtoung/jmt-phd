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
my $sites; ### list of rdd sites (chrom, position, strand, ref_base, rdd_base)
my $tag;
my $alnDB;
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
	"sites=s" => \$sites,
	"tag=s" => \$tag,
	"alnDB=s" => \$alnDB,
	"strand_specific=i" => \$strand_specific,
	"adapter_only=s" => \$adapter_only, 
	"unique_aln_only=i" => \$unique_aln_only,
	"unique_seq_only=s" => \$unique_seq_only,
	"minqual=i" => \$minqual,
	"min_rdd_level=s" => \$min_rdd_level,
	"min_total_count=s" => \$min_total_count
);

################################################################################
### This script gets info for sites that were called RDD elsewhere
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
my $output = $rdd_dir . "/$bam_name.$tag.rdd";
open(OUTPUT,">".$output) or die "[STDERR]: can't open $output: $!\n";
select(OUTPUT);
################################################################################

### LOAD THE LIST OF SITES INTO HASH ###########################################
open(SITES,$sites) or die "[STDERR]: can't open $sites: $!\n";
my %SITES;
while(<SITES>) {
	chomp;
	my ($chrom, $position, $strand, $ref_base, $rdd_base) = split('\t');
	$SITES{$chrom}{$position}{'strand'} = $strand;
	$SITES{$chrom}{$position}{'ref_base'} = $ref_base;
	$SITES{$chrom}{$position}{'rdd_base'} = $rdd_base;
}
################################################################################

### CALL PILEUP ################################################################
my $pileup = "samtools mpileup -f $index -l $sites $bam |";
print STDOUT "pileup command: $pileup\n";
open(PILEUP,$pileup) or die "[STDERR]: can't open $pileup: $!\n";
while(<PILEUP>) {
	chomp;
	my @split = split('\t');
	$DB::single = 1;
	my $CHROM = $split[0];
	my $POSITION = $split[1];
	
	### LOAD INFO FROM SITES FILE
	my ($sites_strand, $sites_ref_base, $sites_rdd_base);
	if (exists $SITES{$CHROM}{$POSITION}) {
		$sites_strand = $SITES{$CHROM}{$POSITION}{'strand'};
		$sites_ref_base = $SITES{$CHROM}{$POSITION}{'ref_base'};
		$sites_rdd_base = $SITES{$CHROM}{$POSITION}{'rdd_base'};
	} else {
		next;
	}

	### PILEUP DATA
	my $RNA_pileup = PileupData->new(\@split,$minqual);

	### DEFINE STRANDS WE'RE INTERESTED IN
	my @STRANDS;
	if ($strand_specific) { 
		if ($sites_strand eq '+') { push(@STRANDS,['+']); }
		elsif ($sites_strand eq '-') { push(@STRANDS,['-']); }
		else { push(@STRANDS,['+','-']); }
	} else { push(@STRANDS,['+','-']); }
	
	my $STATS_OBJ;

	foreach my $STRANDS (@STRANDS) {

		my $TOTAL_COUNT = $RNA_pileup->get_base_count($STRANDS,['A','C','G','T']);
	
		if ($min_total_count) { next if $TOTAL_COUNT < $min_total_count; }

		my $BASE_COUNT = $RNA_pileup->get_base_count($STRANDS,[$sites_rdd_base]);

		if ($min_rdd_level) { next if $TOTAL_COUNT == 0 || $BASE_COUNT/$TOTAL_COUNT < $min_rdd_level; }

		my @RESULTS;

		### CHROM & POSITION
		push(@RESULTS,$CHROM,$POSITION);

		### STRAND
		if (scalar(@{$STRANDS}) == 2) { push(@RESULTS,'*'); }
		elsif ($STRANDS->[0] eq '+') { push(@RESULTS,'+'); }
		else { push(@RESULTS,'-'); }
			
		### REF & RDD BASE
		push(@RESULTS,$sites_ref_base,$sites_rdd_base);

		$STATS_OBJ = make_RDD_stats_object($BAM,$CHROM,$POSITION,$sites_ref_base,$minqual) if scalar(keys %{$STATS_OBJ}) == 0;

		foreach my $ADAPTER_ONLY (split(',',$adapter_only)) {
			foreach my $UNIQUE_SEQ_ONLY (split(',',$unique_seq_only)) {			

				### STATS FOR ALL BASES
				my ($COUNT,$QUAL,$MAPQUAL,$POS5,$POS3) = get_RDD_stats($STATS_OBJ,['A','C','G','T'],$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY);

				### STATS FOR RDD BASE
				my ($RDD_COUNT,$RDD_QUAL,$RDD_MAPQUAL,$RDD_POS5,$RDD_POS3) = get_RDD_stats($STATS_OBJ,[$sites_rdd_base],$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY);

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

sub calc_stats {
	my $ARRAY = shift;

	my $STAT = Statistics::Descriptive::Full->new();
	$STAT->add_data($ARRAY);
	return [$STAT->mean(), $STAT->standard_deviation, $STAT->max(), $STAT->min()];
}
