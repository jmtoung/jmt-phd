#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", 'gpfs/fs121/h/toung/oldhome/dev', '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Database;
use Bio::DB::Sam;
use File::Basename;
use Statistics::Descriptive;
use RegularList;
use PileupData;
use CalculateRDDStats_v2;
use ComplementBase;

my $HOME;
my $sites;
my $region;
my $bam;
my $index;
my $tag; 
my $alnDB;
my $strand_specific = 1;
my $unique_aln_only = 1;
my $unique_seq_only = "0,1";
my $adapter_only = "0";
my $minqual = 20;
my $alnID;

$|++;
my $result = GetOptions(
	"home=s" => \$HOME,
	"sites=s" => \$sites,
	"region=s" => \$region,
	"bam=s" => \$bam,
	"index=s" => \$index,
	"tag=s" => \$tag,
	"alnDB=s" => \$alnDB,
	"strand_specific=i" => \$strand_specific,
	"unique_aln_only=i" => \$unique_aln_only,
	"unique_seq_only=s" => \$unique_seq_only,
	"adapter_only=s" => \$adapter_only, 
	"minqual=i" => \$minqual,
	"alnID=i" => \$alnID
);

################################################################################
### This script extracts RDD info for a bam file given a list of sites & RDD types
################################################################################

### PRINT OPTIONS ##############################################################
print "home:\t$HOME\n";
print "sites:\t$sites\n";
print "region:\t$region\n" if defined $region;
print "bam:\t$bam\n";
print "index:\t$index\n";
defined $tag && print "tag:\t$tag\n" or die "[STDERR]: not defined tag\n"; 
defined $alnDB && print "alnDB:\t$alnDB\n"; # or die "[STDERR]: not defined alnDB\n"; ## changed on 3/31/2012
($strand_specific == 0 || $strand_specific == 1) && print "strand_specific:\t$strand_specific\n" or die "[STDERR]: not defined strand_specific\n";
($unique_aln_only == 0 || $unique_aln_only == 1) && print "unique_aln_only:\t$unique_aln_only\n" or die "[STDERR]: not defined unique_aln_only\n";
($unique_seq_only eq "0,1" || $unique_seq_only eq "0" || $unique_seq_only eq "1") && print "unique_seq_only:\t$unique_seq_only\n" or die "[STDERR]: not defined unique_seq_only\n";
($adapter_only eq "0,1" || $adapter_only eq "0" || $adapter_only eq "1") && print "adapter_only:\t$adapter_only\n" or die "[STDERR]: not defined adapter_only\n";
defined $minqual && print "minqual:\t$minqual\n" or die "[STDERR]: not defined minqual\n";
################################################################################

### GET ALN ID #################################################################
unless (defined $alnID) {
	my $bam_rm_unique = $bam; $bam_rm_unique =~ s/\.unique//; $bam_rm_unique =~ s/\.rmspl//;
	$alnID = Database->new($alnDB)->lookup(1,$bam_rm_unique,0);
	defined $alnID or die "[STDERR]: alnID for '$bam_rm_unique' not defined in '$alnDB'\n";
}
################################################################################

### LOAD BAM & INDEX ###########################################################
my $full_bam; ### added 3/31
if (substr($bam,0,1) eq '/') { $full_bam = $bam; } 
else { $full_bam = $HOME . "/" . $bam; }
print "full_bam:\t$full_bam\n";
-e $full_bam && -e $index or die "[STDERR]: 'bam' $bam or 'index' $index don't exist\n"; 
($full_bam =~ /unique/ or die "[STDERR]: specified unique_aln_only but bam !~ /unique/\n") if $unique_aln_only;
my $BAM = Bio::DB::Sam->new(-bam => $full_bam,-fasta => $index,-autoindex => 1,-split_splices => 1);
################################################################################

### CHECK IF REGION IS DEFINED #################################################
my ($region_index,$region_start,$region_end);
if (defined $region) {
	if ($region =~ /([0-9]+):([0-9]+)-([0-9]+)/) { $region_index = $1; $region_start = $2; $region_end = $3; } 
	else { die "improper region $region\n"; }
}
################################################################################

### OPEN OUTPUT FILE ###########################################################
my ($bam_name,$bam_dir,$bam_exit) = fileparse($full_bam,'\.bam');
$bam_dir =~ s/aln/rdd/;
$bam_name =~ s/\.rmspl//;
my $output = $bam_dir . $bam_name . ".$tag" . "/";
system("mkdir -p $output") unless (-d $output);
$output .= $bam_name . ".$region_index" if defined $region;
$output .= ".rdd";
open(OUTPUT,">$output") or die "[STDERR]: can't open $output: $!\n";
select(OUTPUT);
$|++;
################################################################################

### DETERMINE STRAND SPECIFICITY ###############################################
($strand_specific == 0 or $strand_specific == 1) or die "[STDERR]: define strand specificity\n";
my @strands;
if ($strand_specific) { push(@strands,['+'],['-']); }
else { push(@strands,['+','-']); }
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

	next if $split[0] eq 'chrom';

	my $CHROM = $split[0];
	my $POSITION = $split[1];

	my $STRANDS = [$split[2]];
	$STRANDS = '+,-' if $STRANDS eq '.'; ### HACK SINCE WE'RE USING THE RDD FORMAT NOT NONREF/RDD SITES FORMAT
	my $REF_BASE = $split[3];
	my $RDD_BASE = $split[5];
	$RDD_BASE = $split[4] if ($RDD_BASE eq '1' or $RDD_BASE eq '0'); ### HACK SINCE WE'RE USING RDD FORMAT NOT NONREF/RDD SITES FORMAT
	
	my $RDD_OBJ = make_RDD_stats_object($BAM,$CHROM,$POSITION,$minqual);	

	foreach my $STRANDS (@strands) {
		my ($REF_BASE_new,$RDD_BASE_new);
		if ($STRANDS->[0] eq '-') { $REF_BASE_new = complement_base($REF_BASE); $RDD_BASE_new = complement_base($RDD_BASE); }
		else { $REF_BASE_new = $REF_BASE; $RDD_BASE_new = $RDD_BASE; }
		
		foreach my $ADAPTER_ONLY (split(',',$adapter_only)) {
			foreach my $UNIQUE_SEQ_ONLY (split(',',$unique_seq_only)) {
				my @RESULTS;
				
				### CHROM AND POSITION
				push(@RESULTS,$CHROM,$POSITION);
				
				### STRAND
				if (scalar(@{$STRANDS}) == 2) { push(@RESULTS,'.'); }
				elsif ($STRANDS->[0] eq '+') { push(@RESULTS,'+'); }
				elsif ($STRANDS->[0] eq '-') { push(@RESULTS,'-'); }
				else { die "[STDERR]: invalid strand!!!\n"; }

				### REF AND RDD BASE
				push(@RESULTS,$REF_BASE_new,$RDD_BASE_new);
	
				### GET ALL BASES
				my $TOTAL_BASES = get_bases($RDD_OBJ,$STRANDS);

				### STATS FOR ALL BASES
				my ($COUNT,$QUAL,$MAPQUAL,$POS5,$POS3,$SPLICE,$INSERTION,$DELETION,$MATCH) = get_RDD_stats($RDD_OBJ,$TOTAL_BASES,$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY);
		
				### STATS FOR RDD BASE
				my ($RDD_COUNT,$RDD_QUAL,$RDD_MAPQUAL,$RDD_POS5,$RDD_POS3,$RDD_SPLICE,$RDD_INSERTION,$RDD_DELETION,$RDD_MATCH) = get_RDD_stats($RDD_OBJ,[$RDD_BASE_new],$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY);
				
				### PRINT RESULTS
				print join("\t",@RESULTS), "\t", $unique_aln_only, "\t", $UNIQUE_SEQ_ONLY, "\t", $ADAPTER_ONLY, "\t";
				print $COUNT, "\t", "{";
				print join(",",@{calc_stats($QUAL)}) if scalar(@{$QUAL}) != 0;
				print "}", "\t", "{";
				print join(",",@{calc_stats($MAPQUAL)}) if scalar(@{$MAPQUAL}) != 0;
				print "}", "\t", "{";
				print join(",",@{calc_stats($POS5)}) if scalar(@{$POS5}) != 0;
				print "}", "\t", "{";
				print join(",",@{calc_stats($POS3)}) if scalar(@{$POS3}) != 0;
				print "}", "\t";
				print $SPLICE, "\t", $INSERTION, "\t", $DELETION, "\t", $MATCH, "\t";

				print $RDD_COUNT, "\t", "{";
				print join(",",@{calc_stats($RDD_QUAL)}) if scalar(@{$RDD_QUAL}) != 0;
				print "}", "\t", "{";
				print join(",",@{calc_stats($RDD_MAPQUAL)}) if scalar(@{$RDD_MAPQUAL}) != 0;
				print "}", "\t", "{";
				print join(",",@{calc_stats($RDD_POS5)}) if scalar(@{$RDD_POS5}) != 0;
				print "}", "\t", "{";
				print join(",",@{calc_stats($RDD_POS3)}) if scalar(@{$RDD_POS3}) != 0;
				print "}", "\t";
				print $RDD_SPLICE, "\t", $RDD_INSERTION, "\t", $RDD_DELETION, "\t", $RDD_MATCH, "\t";
					
				print $alnID, "\n";
			}
		}
	}
}
close(SITES);


################################################################################

sub calc_stats {
	my $ARRAY = shift;

	my $STAT = Statistics::Descriptive::Full->new();
	$STAT->add_data($ARRAY);
	return [$STAT->mean(), $STAT->standard_deviation, $STAT->max(), $STAT->min()];
}
