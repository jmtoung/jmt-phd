#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", 'gpfs/fs121/h/toung/oldhome/dev', '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use GetLineCount;
use Database;
use Bio::DB::Sam;
use File::Basename;
use Statistics::Descriptive;
use CalculateRDDStats_v5;
use Data::Dumper;

umask 0007;

my $HOME; ### optional (if bam is not full path)
my $sites; ### REQUIRED
my $region; ### optional
my $bam; ### if you are going to look up alnID, don't give full path (but then you must give HOME)
my $index; ### REQUIRED
my $tag; ### REQUIRED
my $alnDB; ### optional if you are going to look up alnID
my $strand_specific; ### REQUIRED
my $unique_aln_only; ### REQUIRED
my $unique_seq_only; ### REQUIRED
my $adapter_only; ### REQUIRED
my $minqual = 20; ### REQUIRED
my $alnID; ### optional (can put anything; if defined won't look up alnID)
my $combine_trim = 0;
my $stdout = 0;
my $discount_indels = 0;

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
	"alnID=s" => \$alnID,
	"combine_trim=i" => \$combine_trim,
	"stdout=i" => \$stdout,
	"discount_indels=i" => \$discount_indels
);

################################################################################
### This script extracts RDD info for a bam file given a list of sites in bed format
################################################################################

### PRINT OPTIONS ##############################################################
print STDERR "home:\t$HOME\n" if defined $HOME;
(print STDERR "sites:\t$sites\n") && -e $sites or die "[STDERR]: $sites not defined\n";
print STDERR "region:\t$region\n" if defined $region;
(print STDERR "index:\t$index\n") && -e $index or die "[STDERR]: index $index not defined\n";
unless ($stdout) { (print STDERR "tag:\t$tag\n") && defined $tag or die "[STDERR]: not defined tag\n"; }
print STDERR "alnDB:\t$alnDB\n" if defined $alnDB;
(print STDERR "strand_specific:\t$strand_specific\n") && ($strand_specific == 0 || $strand_specific == 1) or die "[STDERR]: not defined strand_specific\n";
(print STDERR "unique_aln_only:\t$unique_aln_only\n") && ($unique_aln_only == 0 || $unique_aln_only == 1) or die "[STDERR]: not defined unique_aln_only\n";
(print STDERR "unique_seq_only:\t$unique_seq_only\n") && ($unique_seq_only eq "0,1" || $unique_seq_only eq "0" || $unique_seq_only eq "1") or die "[STDERR]: not defined unique_seq_only\n";
(print STDERR "adapter_only:\t$adapter_only\n") && ($adapter_only eq "0,1" || $adapter_only eq "0" || $adapter_only eq "1") or die "[STDERR]: not defined adapter_only\n";
(print STDERR "minqual:\t$minqual\n") && $minqual =~ /^[0-9]+$/ or die "[STDERR]: not defined minqual\n";
(print STDERR "combine_trim:\t$combine_trim\n") && ($combine_trim == 0 || $combine_trim == 1) or die "[STDERR]: combine trim not defined\n";
(print STDERR "stdout:\t$stdout\n") && $stdout =~ /^[0-9]+$/ or die "[STDERR]: stdout not defined\n";
(print STDERR "discount_indels:\t$discount_indels\n") && $discount_indels =~ /^[0-9]+$/ or die "[STDERR]: discount indels not defined\n";
################################################################################

### LOAD BAM & INDEX ###########################################################
my @bam = split(',',$bam);
my @BAM; ### bam objects
my @full_bam; ### full bam paths
foreach my $b (@bam) {
	my $full_bam; 
	if (substr($b,0,1) eq '/') { $full_bam = $b; } 
	else { $full_bam = $HOME . "/" . $b; }
	print STDERR "full_bam:\t$full_bam\n";
	($full_bam =~ /unique/ or die "[STDERR]: specified unique_aln_only but bam !~ /unique/\n") if $unique_aln_only;
	push(@full_bam,$full_bam);
	my $BAM = Bio::DB::Sam->new(-bam => $full_bam,-fasta => $index,-autoindex => 1);
	push(@BAM,$BAM);
}
print STDERR "numBamFiles:\t", scalar(@bam), "\n";
################################################################################

### GET ALN ID #################################################################
unless (defined $alnID) {
	my @alnID;
	foreach my $b (@bam) {
		my $bam_rm_unique = $b; $bam_rm_unique =~ s/\.unique//; $bam_rm_unique =~ s/\.rmspl//;
		$alnID = Database->new($alnDB)->lookup(1,$bam_rm_unique,0);
		defined $alnID or die "[STDERR]: alnID for '$bam_rm_unique' not defined in '$alnDB'\n";
		push(@alnID,$alnID);
	}
	$alnID = join(",",sort {$a<=>$b} @alnID);
}
print STDERR "alnID:\t$alnID\n";
################################################################################

### CHECK IF REGION IS DEFINED #################################################
my ($region_index,$region_start,$region_end);
if (defined $region) {
	if ($region =~ /([0-9]+):([0-9]+)-([0-9]+)/) { $region_index = $1; $region_start = $2; $region_end = $3; } 
	else { die "improper region $region\n"; }
}
################################################################################

### OPEN OUTPUT FILE ###########################################################
my $FH = *STDOUT;
unless ($stdout) {
	my ($bam_name,$bam_dir,$bam_ext) = fileparse($full_bam[0],'\.bam');
	$bam_dir =~ s/aln/rdd/;
	$bam_name =~ s/\.rmspl//;
	my $output = $bam_dir . $bam_name . ".$tag" . "/";
	runCommand("mkdkir -p $output") unless -d $output;	
	$output .= $bam_name . "." . 0 x (5 - length($region_index)) . $region_index if defined $region;
	$output .= ".rdd";
	print STDERR "output:\t$output\n";
	open($FH,">$output") or die "[STDERR]: can't open $output: $!\n";
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

	my $RDD_OBJ = makeRddObject(\@BAM,$CHROM,$END,$minqual,$combine_trim);	

	foreach my $STRANDS (@strands) {

		foreach my $ADAPTER_ONLY (split(',',$adapter_only)) {

			foreach my $UNIQUE_SEQ_ONLY (split(',',$unique_seq_only)) {

				my @RESULTS;				
				push(@RESULTS,$CHROM,$START,$END,join(",",@{$STRANDS}),$split[4],$split[5]);
				push(@RESULTS,$unique_aln_only,$UNIQUE_SEQ_ONLY,$ADAPTER_ONLY);

				### GET ALL BASES
				my $TOTAL_BASES = getBases($RDD_OBJ,$STRANDS,$discount_indels);

				### DO THE FOLLOWING FOR TWO GROUPS OF BASES
				foreach my $GROUP ($TOTAL_BASES,[$RDD_BASE]) {

					my $STATS = getRddStats($RDD_OBJ,$GROUP,$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY);
					foreach my $STAT (@{$STATS}) {
						if (ref($STAT) eq '') {
							push(@RESULTS,$STAT);
						} elsif (ref($STAT) eq 'ARRAY') {
							push(@RESULTS,join(";",@{calcStats($STAT)}));
						} else {
							die "[STDERR]: weird referenced object!\n";
						}
					}
				}
				
				push(@RESULTS,$alnID);
				
				print {$FH} join("\t",@RESULTS), "\n";
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
		return [$STAT->mean(), $STAT->standard_deviation, $STAT->max(), $STAT->min()];
	} else {
		return ['NA','NA','NA','NA'];
	}
}
