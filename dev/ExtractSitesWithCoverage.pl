#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Database;
use File::Basename;
use PileupData;
use IlmnDNA;
use ComplementBase;
umask 0007;

my $home;
my $pileup;
my $region; ### added on 4/17/2012
my $strand_specific;
my $minqual;
my $min_total_reads = 10;

#############################################################################################

my $result = GetOptions(
	"home=s" => \$home,
	"pileup=s" => \$pileup,
	"region=s" => \$region,
	"strand_specific=i" => \$strand_specific,
	"minqual=i" => \$minqual,
	"min_total_reads=i" => \$min_total_reads
);

################################################################################
### This script outputs sites where RNA is different from DNA subject to minqual thresholds
### the 'bed' option means to just spit out the sites that were not skipped so we can tally up how many gigabases of sequenced was analyzed
################################################################################

(print STDERR "home:\t$home\n") && -d $home or die "[STDERR]: home not defined\n";
(print STDERR "pileup:\t$pileup\n") && -e $pileup or die "[STDERR]: pileup not defined\n";
(print STDERR "region:\t$region\n") if defined $region;
(print STDERR "strand_specific:\t$strand_specific\n") && defined $strand_specific or die "[STDERR]: strand specific not defined\n";
(print STDERR "minqual:\t$minqual\n") && defined $minqual or die "[STDERR]: minqual not defined\n";
(print STDERR "min_total_reads:\t$min_total_reads\n") or die "[STDERR]: min_total_reads not defined\n";

my ($regionPart,$regionStart,$regionEnd);
if (defined $region) {
	($regionPart,$regionStart,$regionEnd) = split('-',$region);
}

my ($pileup_name,$pileup_dir,$pileup_ext) = fileparse($pileup,'\.pileup(\.gz)?');
my $PILEUP;
if (defined $regionStart) {
	$PILEUP = PileupData->new($pileup,$minqual,$regionStart);
} else {
	$PILEUP = PileupData->new($pileup,$minqual);
}

my @strands;
if ($strand_specific) { push(@strands,['+'],['-']); }
else { push(@strands,['+','-']); }

while($PILEUP->has_next) {
	my $chrom = $PILEUP->get_chrom; 
	my $position = $PILEUP->get_position;

	if (defined $region) {
		last if $position > $regionEnd;
	}
	
	my $ref_base = $PILEUP->get_ref_base;

	goto NEXT if $ref_base eq 'N';


	foreach my $strands (@strands) {
		my $bases = $PILEUP->get_bases($strands);
		
		my @total_bases;
		foreach my $base (@{$bases}) {
			push(@total_bases,$base) unless $base eq 'S';
		}
		my $total_count = $PILEUP->get_base_count($strands,\@total_bases);

		if ($total_count >= $min_total_reads) {
			print join("\t",$chrom, $position - 1, $position, join(",",@{$strands})), "\n";
			goto NEXT;
		}
	}
	
	NEXT: 

	$PILEUP->parse_next_line;
}
