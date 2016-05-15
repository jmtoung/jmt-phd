#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use File::Basename;
use RegularList;
use PileupData;
use ComplementBase;
use List::Util qw[min max];

my $pileup = "/home/jmtoung/Lab/aln/twins/GM14381_s_7_sequence_trimlowqual/bowtie.GM14381_s_7_sequence_trimlowqual.n2e120.unique.pileup/bowtie.GM14381_s_7_sequence_trimlowqual.n2e120.unique.1.pileup,/home/jmtoung/Lab/aln/twins/GM14382_s_8_sequence_trimlowqual/bowtie.GM14382_s_8_sequence_trimlowqual.n2e120.unique.pileup/bowtie.GM14382_s_8_sequence_trimlowqual.n2e120.unique.1.pileup";
my $tag;
my $strand_specific; ### 1 (yes) or 0 (no)
my $minqual;

my $result = GetOptions(
	"pileup=s" => \$pileup,
	"tag=s" => \$tag,
	"strand_specific=i" => \$strand_specific,
	"minqual=i" => \$minqual
);

### LOAD PILEUP FILES ##########################################################
my %PILEUP;
foreach my $P (split(',',$pileup)) { 	
	-e $P or print STDERR "[STDERR]: $P doesn't exist\n";
	$PILEUP{$P} = RegularList->new($P) if -e $P;
}
################################################################################

################################################################################
while(has_next(\%PILEUP)) {
	
	my $min_P = get_min_P(\%PILEUP);

	foreach my $min_P (@{$min_P}) {
		
		my $P_line = $PILEUP{$min_P}->get_next_line;
		my $parsed_P = PileupData->new($P_line,$minqual);
		$DB::single = 1;
		print ".";
		
		
	}
#	my %bases;
#	for(my $i=0; $i <= $#PILEUP; $i++) {
#		if ($PILEUP[$i]->get_column(1) == $PILEUP[$min]->get_column(1)) {
#			my $RNA_pileup = PileupData->new(\@split,$minqual);

#			%bases	
#			$PILEUP[$i]->update_next_line;
#		}
#	}
	
	### print out stuff
	
	
	
}

sub has_next {
	my $PILEUP = shift;
	
	foreach my $P (keys %{$PILEUP}) {
		return 1 if $PILEUP->{$P}->has_next;
	}
	return 0;
}

sub get_min_P {
	my $PILEUP = shift;
	
	my @positions;
	foreach my $P (keys %{$PILEUP}) {
		next unless $PILEUP->{$P}->has_next;
		my $position = $PILEUP->{$P}->get_column(1);
		push(@positions,$position);	
	}
	
	my $min_position = min(@positions);
	
	my @P;
	foreach my $P (keys %{$PILEUP}) {
		push(@P,$P) if ($PILEUP->{$P}->get_column(1) == $min_position);
	}
	
	return \@P;
}

#	chomp;
#	my @split = split('\t');

#	my $CHROM = $split[0];
#	my $POSITION = $split[1];

#	next unless ($CHROM eq $region_chrom && $POSITION >= $region_start && $POSITION <= $region_end);

#	foreach my $EXCLUDE (@EXCLUDE) {
#		while($EXCLUDE->has_next && $EXCLUDE->get_column(3) < $POSITION) { $EXCLUDE->update_next_line; }
#		goto NEXT if ($EXCLUDE->has_next && $EXCLUDE->get_column(2) <= $POSITION);
#	}
#	
#	my $RNA_pileup = PileupData->new(\@split,$minqual);
#	goto NEXT unless $RNA_pileup->is_def;

#	my @STRANDS;
#	if ($strand_specific) { push(@STRANDS,['+'],['-']); }
#	else { push(@STRANDS,['+','-']); }

#	my $REF_BASE = $split[2];
#	my $STATS_OBJ;

#	foreach my $STRANDS (@STRANDS) {

#		my $TOTAL_COUNT = $RNA_pileup->get_base_count($STRANDS,['A','C','G','T']);
#		next if $TOTAL_COUNT == 0;
#		if ($min_total_count) { next if $TOTAL_COUNT < $min_total_count; }

#		foreach my $BASE ('A','C','T','G') {
#			my @RESULTS;

#			if ($STRANDS->[0] eq '-') { next if $BASE eq complement_base($REF_BASE); }
#			else { next if $BASE eq $REF_BASE; }

#			my $BASE_COUNT = $RNA_pileup->get_base_count($STRANDS,[$BASE]);
#			next if $BASE_COUNT == 0;
#			if ($min_rdd_level) { next if $BASE_COUNT/$TOTAL_COUNT < $min_rdd_level; }

#			### CHROM & POSITION
#			push(@RESULTS,$CHROM,$POSITION);

#			### STRAND
#			if (scalar(@{$STRANDS}) == 2) { push(@RESULTS,'*'); }
#			elsif ($STRANDS->[0] eq '+') { push(@RESULTS,'+'); }
#			else { push(@RESULTS,'-'); }
#			
#			### REF & RDD BASE
#			if ($STRANDS->[0] eq '-') { push(@RESULTS,complement_base($REF_BASE),$BASE); }
#			else { push(@RESULTS,$REF_BASE,$BASE); }

#			$STATS_OBJ = make_RDD_stats_object($BAM,$CHROM,$POSITION,$REF_BASE,$minqual) if scalar(keys %{$STATS_OBJ}) == 0;

#			foreach my $ADAPTER_ONLY (split(',',$adapter_only)) {
#				foreach my $UNIQUE_SEQ_ONLY (split(',',$unique_seq_only)) {			

#					### STATS FOR ALL BASES
#					my ($COUNT,$QUAL,$MAPQUAL,$POS5,$POS3) = get_RDD_stats($STATS_OBJ,['A','C','G','T'],$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY);

#					### STATS FOR RDD BASE
#					my ($RDD_COUNT,$RDD_QUAL,$RDD_MAPQUAL,$RDD_POS5,$RDD_POS3) = get_RDD_stats($STATS_OBJ,[$BASE],$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY);

#					next if $COUNT == 0 && $RDD_COUNT == 0;

#					### PRINT RESULTS
#					print join("\t",@RESULTS), "\t", $unique_aln_only, "\t", $UNIQUE_SEQ_ONLY, "\t", $ADAPTER_ONLY, "\t";

#					print $COUNT, "\t";
#					print "{";
#					print join(",",@{calc_stats($QUAL)}) if scalar(@{$QUAL}) != 0;
#					print "}", "\t";
#					print "{";
#					print join(",",@{calc_stats($MAPQUAL)}) if scalar(@{$MAPQUAL}) != 0;
#					print "}", "\t";
#					print "{";
#					print join(",",@{calc_stats($POS5)}) if scalar(@{$POS5}) != 0;
#					print "}", "\t";
#					print "{";
#					print join(",",@{calc_stats($POS3)}) if scalar(@{$POS3}) != 0;
#					print "}", "\t";

#					print $RDD_COUNT, "\t";
#					print "{";
#					print join(",",@{calc_stats($RDD_QUAL)}) if scalar(@{$RDD_QUAL}) != 0;
#					print "}", "\t";
#					print "{";
#					print join(",",@{calc_stats($RDD_MAPQUAL)}) if scalar(@{$RDD_MAPQUAL}) != 0;
#					print "}", "\t";
#					print "{";
#					print join(",",@{calc_stats($RDD_POS5)}) if scalar(@{$RDD_POS5}) != 0;
#					print "}", "\t";
#					print "{";
#					print join(",",@{calc_stats($RDD_POS3)}) if scalar(@{$RDD_POS3}) != 0;
#					print "}", "\t";
#					
#					print $bamID, "\n";
#				}
#			}
#		}
#	}
#	NEXT:
##}

#sub calc_stats {
#	my $ARRAY = shift;

#	my $STAT = Statistics::Descriptive::Full->new();
#	$STAT->add_data($ARRAY);
#	return [$STAT->mean(), $STAT->standard_deviation, $STAT->max(), $STAT->min()];
#}
