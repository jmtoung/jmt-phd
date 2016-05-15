#!/usr/bin/perl -w

## parse the output of BLAT genomic sequences around the RDD
## only look at substitutions, no indels
## each RNA variant gets its own fasta entry,

## query name contains the following information (seperated by ";"):
## 1) chrom 2) position 3) flanking distance 4)ref base 5) rdd base 6) type of molecule (RNA or DNA)

## filter blat output:
## criterion:
## 1) the number of mismatches (M) allowed. M = (querySize + 2)/kmer 
## 2) score of alignment (#matches - #mismatches - #gaps in query - #gaps in target) has to be greater than querySize - M
## 3) user specified maximum number of blocks (sub-alignment with no gaps)
## 

use lib '/Users/zhaoyue/projects/adar/src/jonathan/perl_libraries', '/Users/zhaoy/projects/adar/src/jonathan/perl_libraries', '/ifs/h/zhaoyue/projects/adar/src/jonathan/perl_libraries';
use strict;
use Getopt::Long;
use File::Basename;
use POSIX;
use ComplementBase;
my $pslxFile;
my $inputFile;

# defaults to only look at ungapped alignments
my $maxBlocks = 1;
my $debug = 0;
my $gapPen = 1;
my $options = GetOptions(
    "inputFile=s" => \$inputFile,
    "pslxFile=s" => \$pslxFile,
    "maxBlocks=i" => \$maxBlocks,
    "gapPen=i" => \$gapPen,    
    "debug=i" => \$debug
);

-e $pslxFile or die"[STDERR]: $pslxFile doesn't exist\n";
## print "pslxFile = $pslxFile\n";

my ($fileName, $fileDir, $fileExt) = fileparse($pslxFile, "\.pslx");
## print "name: $fileName, dir: $fileDir, ext: $fileExt\n";

## filtered blat hits
my %outputs;
my %blatResults;
## blat score, #matches - #mismatches - #gaps in query - #gaps in target
my $score;
my $pslxLine;
my $kmer = 12; #size of kmer used for blat
my $maxScore;

## variables representing all the information in a line of blat output (pslx, includes matching sequences)
my ($matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts,$qSeq,$tSeq);

## variables representing information stored in query sequence name
my ($chrom,$position,$flankDist,$refBase,$rddBase,$mol) ;

## position of RDD on the query sequence
my $qRDDPos;

## position of RDD in the alignment block
my $blockRDDPos;
 
my $regionStart;

my @blocksLen;
my @qBlockStarts;
my @tBlockStarts;
my @tBlockSeqs;
my $qBlockStart;
my $qBlockEnd;
my $rddBlock;
my $tRDDBase;

my $rddExplained;
my $numMismatches;
my $mismatchAllowed;

my $tSpan;
my $numMM;

## go through pslx file, put good BLAT hits in hash, key is query name, value is array of all blat output
open(FILE, $pslxFile) or die "[STDERR]: can't open $pslxFile:$!\n";
while(<FILE>) {
    chomp;
    $pslxLine = $_;
    ($matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts,$qSeq,$tSeq) = split('\t', $pslxLine);

    ## these fields end with trailing "," chop it off
    chop($qSeq);
    chop($tSeq);
    chop($blockSizes);	
    chop($qStarts);
    chop($tStarts);

    $mismatchAllowed = floor(($qSize + 2)/$kmer);
    # skip if there are any N's in the sequence
#    next if $nCount > 0;

    # skip if there are gapps in thgape query sequence
    # gaps in target is fine, result of aligning a read that spans exon-exon junctions
    # gaps in query is a problem, ignore for now
    #next if $qNumInsert > 0;

    # quick'n dirty filter, get rid of blat hits whose scores are less than half of query length
    # $score = $matches - $misMatches - $qNumInsert - $tNumInsert;
    # next if ($score < ($qSize/2));

    ($chrom,$position,$flankDist,$refBase,$rddBase,$mol) = split(';',$qName);


    ## start of genomic region we used to blat, (subtract 1 to make it 0-based)
    $regionStart = $position - $flankDist -1;
    
    ## ignore hit to the region containing RDD we used to blat in the first place
    ## tStart is start of alignment in the genome
    ## qStart is the start of alignment on the read, alignment may not start from the first base of query 
    next if ($chrom eq $tName && (($tStart - $qStart) == $regionStart));

    ## allow insertions in the target (gapped alignment to genome), but not insertion in query
    next if ($qNumInsert > 0);

    ## 0-based RDD position on the read, forward strand, the same as $flankDist
    $qRDDPos = $flankDist;

    ## qStart and qEnd coordinates are always with regard to forward strand, regardless of which stand query was mapped to
    ## since qRDDPos is also with regard to forward strand, no adjustments are necessary

    ## ignore alignments that ends before RDD position
#    next if $qEnd < $qRDDPos;

    ## ignore alignments that starts after RDD position
#    next if $qStart > $qRDDPos;

    ## ignore alignments with more than max number of blocks
    next if $blockCount > $maxBlocks;

    ## figure out which block contains the RDD site
    @blocksLen = split(",", $blockSizes);
    @qBlockStarts = split(",", $qStarts);
    @tBlockStarts = split(",", $tStarts);
    @tBlockSeqs = split(",", $tSeq);

    $rddBlock = -1;
    for (my $i = 0; $i < $blockCount; $i++) {
	## blockStarts coordinates are with regard to reverse strand if query is mapped to reverse strand
	$qBlockStart = $qBlockStarts[$i];
	$qBlockEnd = $qBlockStart + $blocksLen[$i] - 1;
	if ($strand eq "-") {
	    ## qBlockStart and qBlockEnd are with regard to the reverse strand, make qRDDPos the same way
	    $qRDDPos = $qSize - 1 - $qRDDPos;
	}

	print "$pslxLine\n" if ($debug);
	print "qBlockStart:$qBlockStart\tqBlockEnd:$qBlockEnd\tqRDDPosition:$qRDDPos\n" if ($debug);
	## if RDD position is located within the current block, check to see if the aligned base in target is the same as RNA allele
	if ($qRDDPos >= $qBlockStart && $qRDDPos <= $qBlockEnd) {
	    $rddBlock = $i;
	    $blockRDDPos = $qRDDPos - $qBlockStarts[$i];
	    $tRDDBase = substr($tBlockSeqs[$i], $blockRDDPos, 1);
	    $tRDDBase = complement_base($tRDDBase) if ($strand eq "-");
	    $tRDDBase = uc($tRDDBase);
	    
	    print "blockRDDPos:$blockRDDPos\tRDD Base:$rddBase\ttRDDBase:$tRDDBase\n" if ($debug);

	    ## the current hit explains the RDD (has the same base as RNA allele in the aligned position)
	    if ($tRDDBase eq $rddBase) {
		## score of the hit
#		$score = $matches + $repMatches - $misMatches - $qNumInsert - $tNumInsert;		
		## since we BLAT the RNA form, +1 to account for that mismatch
		$numMM = ($qSize - ($matches + $repMatches));
		$score =  $numMM + $gapPen*$tNumInsert;
		$maxScore = floor(($qSize+2)/$kmer);

		$tSpan = $tEnd - $tStart;

		## hit explains RDD and meets the minimum threashold		
		if ($score <= $maxScore) {
		    $blatResults{$chrom}{$position}{$refBase}{$rddBase}{$flankDist}{$tNumInsert}{$numMM}++;
#		    print "$chrom\t$position\t$tName\t$strand\t$tStart\t$refBase\t$rddBase\t$tNumInsert\t$tRDDBase\t$numMM\t$score\t$maxScore\n"; 
		}
	    } else {
	    }
	    next;
	}
    }
}
close(FILE);

print "Finished Parsing BLAT output\n";

my @numMismatches;
my $out;

open(INFILE, $inputFile) or die "[STDERR]: can't open $inputFile: $!\n";
my ($file_name,$file_dir,$file_ext) = fileparse($inputFile,'\.txt');
my $output = $file_dir . $file_name . "_blat-explRDD_gapped" . $file_ext;
open(OUTFILE,">".$output) or die "[STDERR]: can't open $output: $!\n";

while (<INFILE>) {
    chomp;
    my @split = split('\t');
    my $chrom = $split[0];
    my $position = $split[1];
    my $strand = $split[2];
    die "[STDERR]: strand is not +,-\n" unless $strand eq "+,-";

    my $refBase = $split[3];
    my $rddBase = $split[5];


    $out = "$chrom\t$position\t$refBase\t$rddBase\t";
    for $tNumInsert (0,1,2,3) {
	if (exists $blatResults{$chrom}{$position}{$refBase}{$rddBase}{$flankDist}{$tNumInsert}) {
	    for $numMM (sort {$a<=>$b} keys %{ $blatResults{$chrom}{$position}{$refBase}{$rddBase}{$flankDist}{$tNumInsert} }) {
		if (exists $blatResults{$chrom}{$position}{$refBase}{$rddBase}{$flankDist}{$tNumInsert}{$numMM}) {
		    push @numMismatches, $numMM . ":" . $blatResults{$chrom}{$position}{$refBase}{$rddBase}{$flankDist}{$tNumInsert}{$numMM};
		}
	    }
	    push (@numMismatches, "-") unless (@numMismatches);
	    $out = $out . join(";", @numMismatches) . "\t";
	    
	} else {
	    $out = $out . "-\t";
	}
	@numMismatches = ();
    }
#		    $out = chomp($out);
    print OUTFILE "$out\n";
}
 
close(INFILE);
close(OUTFILE);   

# for $chrom ( keys %blatResults ) {
#     for $position ( keys %{ $blatResults{$chrom} } ) {
# 	for $refBase (keys %{ $blatResults{$chrom}{$position} }) {
# 	    for $rddBase (keys %{ $blatResults{$chrom}{$position}{$refBase} }) {
# 		for $flankDist (keys  %{ $blatResults{$chrom}{$position}{$refBase}{$rddBase} }) {
# 		    $out = "$chrom\t$position\t$refBase\t$rddBase\t";

# 		    for $tNumInsert (0,1,2,3) {
# 			if (exists $blatResults{$chrom}{$position}{$refBase}{$rddBase}{$flankDist}{$tNumInsert}) {
# 			    for $numMM (sort {$a<=>$b} keys %{ $blatResults{$chrom}{$position}{$refBase}{$rddBase}{$flankDist}{$tNumInsert} }) {
# 				if (exists $blatResults{$chrom}{$position}{$refBase}{$rddBase}{$flankDist}{$tNumInsert}{$numMM}) {
# 				    push @numMismatches, $numMM . ":" . $blatResults{$chrom}{$position}{$refBase}{$rddBase}{$flankDist}{$tNumInsert}{$numMM};
# 				}
# 			    }
# 			    push (@numMismatches, "-") unless (@numMismatches);
# 			    $out = $out . join(";", @numMismatches) . "\t";

# 			} else {
# 			    $out = $out . "-\t";
# 			}
# 			@numMismatches = ();
# 		    }
# #		    $out = chomp($out);
# 		    print "$out\n";
# 		}
# 	    }
# 	}
#     }
# }
# }
