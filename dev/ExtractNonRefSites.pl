#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Database;
use File::Basename;
use PileupData;
use RegularList;
use ComplementBase;

my $pileup;
my $tag;
my $exclude;
my $strand_specific;;
my $minqual;
my $alnDB;
my $alnID;

my $result = GetOptions(
	"pileup=s" => \$pileup,
	"tag=s" => \$tag,
	"exclude=s" => \$exclude,
	"strand_specific=i" => \$strand_specific,
	"minqual=i" => \$minqual,
	"alnDB=s" => \$alnDB,
	"alnID=s" => \$alnID
);

################################################################################
### This script calls outputs sites where nonref alleles exist
################################################################################

-e $pileup or die "[STDERR]: can't open $pileup: $!\n";
my ($pileup_name,$pileup_dir,$pileup_ext) = fileparse($pileup,'\.pileup(\.gz)?');

my $output_file = $pileup_dir . $pileup_name . ".$tag.txt";
open(OUTPUT,">$output_file") or die "[STDERR]: can't open $output_file: $!\n";
my $old_fh = select(OUTPUT);
$|++;
select($old_fh);

my $bam = $pileup_dir;
$bam =~ s/(\.unique)?\.pileup\//\.bam/;

### LOOK UP alnID if it doesn't already exit ###################################
unless (defined $alnID) {
	$alnID = Database->new($alnDB)->lookup(1,$bam,0);
	defined $alnID or die "[STDERR]: alnID for '$bam' not defined in '$alnDB'\n";
}
################################################################################

### LOAD sites we want to exclude ##############################################
if ($exclude) { -e $exclude or die "[STDERR]: can't open $exclude: $!\n"; }
my $EXCLUDE = RegularList->new($exclude);
################################################################################

my @strands;
if ($strand_specific) { push(@strands,['+'],['-']); }
else { push(@strands,['+','-']); }

my (%exclude_previous_sources,$exclude_previous_position);
$exclude_previous_position = -1;

################################################################################
my $PILEUP = PileupData->new($pileup,$minqual);
while($PILEUP->has_next) {

	my $chrom = $PILEUP->get_chrom; 
	my $position = $PILEUP->get_position;
	my $ref_base = $PILEUP->get_ref_base;
	
	goto NEXT if $ref_base eq 'N';

	while($EXCLUDE->has_next && $EXCLUDE->get_column(3) < $position) { 

		%exclude_previous_sources = ();
		foreach my $source (split(',',$EXCLUDE->get_column(0))) { $exclude_previous_sources{$source}++; }
		$exclude_previous_position = $EXCLUDE->get_column(2);

		$EXCLUDE->update_next_line;
	}

	### determine if there's a neighboring polymorphism
	my @neighboring_poly;
	if ($EXCLUDE->has_next && $EXCLUDE->get_column(2) <= $position) {
		goto NEXT;
	} else {
		if (($position - $exclude_previous_position) <= 3 && ($position - $exclude_previous_position) > 0) {
			foreach my $source ('1000G-INDEL-Sites','1000G-Pileup-del','1000G-Pileup-ins','snp130-del','snp130-indel','snp130-ins') {
				push(@neighboring_poly,$source) if exists $exclude_previous_sources{$source};
			}
		}
	}
	my $neighboring_poly;
	if (@neighboring_poly) { $neighboring_poly = join(',',@neighboring_poly); }
	else { $neighboring_poly = '-'; }
	###

	foreach my $strands (@strands) {
		my $bases = $PILEUP->get_bases($strands);
		
		my @total_bases;
		foreach my $base (@{$bases}) {
			push(@total_bases,$base) unless $base eq 'S';
		}
		my $total_count = $PILEUP->get_base_count($strands,\@total_bases);

		foreach my $base (@{$bases}) {
			next if $base eq 'S';
			
			my $REF_BASE;
			if ($strands->[0] eq '-') { $REF_BASE = complement_base($ref_base); }
			else { $REF_BASE = $ref_base; }

			next if $base eq $REF_BASE;
			
			my $base_count = $PILEUP->get_base_count($strands,[$base]);

			print OUTPUT $chrom, "\t", $position, "\t", join(",",@{$strands}), "\t", $REF_BASE, "\t", $total_count, "\t", $base, "\t", $base_count, "\t{", $neighboring_poly, "}\t", $alnID, "\n";
		}
	}
	
	NEXT: 
	
	$PILEUP->parse_next_line;
}
close(OUTPUT);
