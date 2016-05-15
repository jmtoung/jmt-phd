#!/usr/bin/perl -w

use strict;
use POSIX;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev', '/gpfs/fs121/h/toung/dev';
use File::Basename;

my $pslxFile;
my $kmer = 12;
my $misMatchPen = 0; ### number of mismatches allowed relative to Gsnap amount (readlength+2)/kmer - 2
my $maxBlocks = 3;
my $gapPen = 0;

my $options = GetOptions(
	"pslxFile=s" => \$pslxFile,
	"kmer=i" => \$kmer,
	"misMatchPen=i" => \$misMatchPen,
	"maxBlocks=i" => \$maxBlocks,
	"gapPen=i" => \$gapPen
);

if (defined $pslxFile) { 
	(print STDERR "pslxFile:\t$pslxFile\n") && -e $pslxFile or die "[STDERR]: pslxFile not defined\nexists\n";
}
(print STDERR "kmer:\t$kmer\n") && $kmer =~ /^[0-9]+$/ or die "[STDERR]: kmer not defined properly\n";
(print STDERR "misMatchPen:\t$misMatchPen\trelative to gsnap\n") && $misMatchPen =~ /^[0-9]+$/ or die "[STDERR]: misMatchPen not defined properly\n";
(print STDERR "maxBlocks:\t$maxBlocks\n") && $maxBlocks =~ /^[0-9]+$/ or die "[STDERR]: undefined maxBlocks\n";
(print STDERR "gapPen:\t$gapPen\n") && $gapPen =~ /^[0-9]+$/ or die "[STDERR]: undefined gapPen\n";

my $FH = *STDIN;
if (defined $pslxFile) { 
	open($FH,$pslxFile) or die "[STDERR]: can't open $pslxFile: $!\n";
}

my %stats;
while(<$FH>) {
	chomp;
	
	my ($matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts,$qSeq,$tSeq) = split('\t');
	$blockSizes =~ s/,$//; $qStarts =~ s/,$//; $tStarts =~ s/,$//; $qSeq =~ s/,$//; $tSeq =~ s/,$//;

	### query info ### note that ref and rdd base are both relative to the strand of the rdd ($rddStrand)
	my ($chrom,$position,$rddStrand,$flank_dist,$ref_base,$rdd_base,$rddId) = split(';',$qName);
	my $region_start = $position - $flank_dist - 1; ### make it 0-based so can compare to $tStart!  ### the start of the region we blatted

	$stats{'0_numRddIds'}{$flank_dist}{$rddId}++;
	$stats{'1_numBlatAln'}{$flank_dist}++;

	### [1] first check if the result is the same as where we yoinked it out of. for this to be true:
	### strand must be on positive strand b/c we took flanking sequence from positive strand
	### chrom ($qName) must equal $tName
	### $region_start (0-based start of region we got it from) must equal the $tStart (minus $qStart if alignment doesn't start at position 0)
	if ($strand eq '+' && $chrom eq $tName && (($tStart - $qStart) == $region_start)) {
		$stats{'2_aln2SameLocation'}{$flank_dist}++;
		next;
	}

	### [2] calculate the number of mismatches (M) allowed. M = (readlength + 2)/kmer -2 + $misMatchPen ($misMatchPen is to allow some wiggle room)
	$qSize == ($flank_dist*2 + 1) or die "[STDERR]: unexpected qSize ($qSize) for flank_dist ($flank_dist)\n";
	my $misMatches_allowed = floor(($qSize + 2)/$kmer - 2) + $misMatchPen;

	my $numMisMatches = $qSize - $matches - $repMatches + $gapPen * ($blockCount - 1); ### cannot use $misMatches b/c it doesn't count things not aligned as mismatches

	if ($numMisMatches > $misMatches_allowed) {
		$stats{'3_tooManyMismatches'}{$flank_dist}++;
		next;
	}
		
	### [3] require block count be no greater than $maxBlocks
	if ($blockCount > $maxBlocks) {
		$stats{'4_tooManyBlockCount'}{$flank_dist}++;
		next;
	}

	### [4] find out if RDD position/base is in alignment AND whether it matches target or not
	my $qRddPos = $flank_dist; ### 0-based distance of RDD position in query

	## qBlockStart ($qStarts) and blockSizes are reerse ordered, make qRDDPos the same way
	$qRddPos = $qSize - 1 - $qRddPos if $strand eq '-';

	my @blockSizes = split(',',$blockSizes);
	my @qStarts = split(',',$qStarts);
	my @tStarts = split(',',$tStarts);
	my @qSeq = split(',',$qSeq);
	my @tSeq = split(',',$tSeq);
	
	check_equal_size([\@blockSizes,\@qStarts,\@tStarts,\@qSeq,\@tSeq]); ### check that the arrays are all the same size

	for (my $i = 0; $i < @blockSizes; $i++) {
		my $qStart = $qStarts[$i];
		my $qEnd = $qStart + $blockSizes[$i] - 1;
	
		if ($qRddPos >= $qStart && $qRddPos <= $qEnd) {
			my $blockRddPos = $qRddPos - $qStarts[$i]; ### position of Rdd in current block

			### the base in the Rdd position in the target
			my $tRddBase = substr($tSeq[$i],$blockRddPos,1);
#######################	$tRddBase =~ tr/ACGT/TGCA/ if ($strand eq '-'); ### not necessary since we'll do the same for $qRddBase below
			$tRddBase = uc($tRddBase);

			### the base in the Rdd position in the query
			my $qRddBase = substr($qSeq[$i],$blockRddPos,1);
#######################	$qRddBase =~ tr/ACGT/TGCA/ if ($strand eq '-'); ### not necessary, see above
			$qRddBase = uc($qRddBase);
			
			my $rddExplained = 0;			
			$rddExplained = 1 if $tRddBase eq $qRddBase;

			print join("\t",$rddId,$flank_dist,$rddExplained,$numMisMatches,$misMatches_allowed,$tName,$tStarts,$blockSizes,$strand,$chrom,$rddStrand,$position,$ref_base,$rdd_base), "\n"; 

			$stats{'5_numRddContainedInAln'}{$flank_dist}++;
		}
	}
					
}	

close($FH);

$pslxFile = '-' if !defined $pslxFile;

foreach my $key (sort keys %stats) {
	foreach my $flank_dist (sort {$a<=>$b} keys %{$stats{$key}}) {
		if ($key eq '0_numRddIds') {
			my @num = keys %{$stats{$key}{$flank_dist}};
			print STDERR join("\t",$key,$flank_dist,scalar(@num),$pslxFile), "\n";
		} else {
			print STDERR join("\t",$key,$flank_dist,$stats{$key}{$flank_dist},$pslxFile), "\n";
		}
	}
}

sub check_equal_size {
	my $arrays = shift;
	
	my %sizes;
	foreach my $array (@{$arrays}) {
		$sizes{scalar(@{$array})}++;
	}
	my @sizes = keys %sizes;
	scalar(@sizes) == 1 or die "[STDERR]: array sizes not equal\n";
}
