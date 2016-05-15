#!/usr/bin/perl -w

use lib '/gpfs/fs121/h/toung/oldhome/dev';
use Database;
use strict;
use Getopt::Long;

my $files;
my $alnDB;
my $fastqDB;
my $samplesDB;

my $options = GetOptions(
	"files=s" => \$files,
	"alnDB=s" => \$alnDB,
	"fastqDB=s" => \$fastqDB,
	"samplesDB=s" => \$samplesDB
);

my @files = split(',',$files); # = readdir(DIR);

my %results;
### relation stores whether a comparison is a twin or unrelated
### other stores the alnID of the other sample
my (%relation,%other);

foreach my $file (@files) {
	
	open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";

	while(<FILE>) {
		chomp;
		my ($comparison,$chrom,$num_samples,$sample,$num_genotypes,$genotypes,$zygosity,$count) = split('\t');
		my $type = calc_type($genotypes);
	
		### do this once only...calculate %relation and %other for this $comparison
		unless (defined $relation{$comparison}) {
			### look up the family IDs of each person in the comparison
			my %family;
			my @alnID = split(',',$comparison); ### the alnID in the comparison
			### define the other alnIDs for each alnID in this comparison
			for (my $i=0; $i < @alnID; $i++) {
				for (my $j=0; $j < @alnID; $j++) {
					next if $i==$j; 
					$other{$comparison}{$alnID[$i]}{$alnID[$j]}++;
				}
			}
			foreach my $alnID (sort {$a<=>$b} keys %{$other{$comparison}}) {
				$other{$comparison}{$alnID} = join(",",sort {$a<=>$b} keys %{$other{$comparison}{$alnID}});
			}
			### determine the relation for this comparison
			foreach my $alnID (@alnID) {
				my $fastqID = Database->new($alnDB)->lookup(0,$alnID,2);
				defined $fastqID or die "[STDERR]: can't find fastqID for $alnID in alnDB '$alnDB'\n";
				my $sample = Database->new($fastqDB)->lookup(0,$fastqID,3);
				defined $sample or die "[STDERR]: can't find sample for $fastqID in fastqDB '$fastqDB'\n";
				my $family = Database->new($samplesDB)->lookup(0,$sample,2);
				defined $family or die "[STDERR]: can't fine family for $sample in samplesDB '$samplesDB'\n";
				$family{$family}++;
			}
			if (keys %family == 1) {
				$relation{$comparison} = 'twins';
			} else {
				$relation{$comparison} = 'unrelateds';
			}
		}

		### add to expression concordance 
		foreach my $S (split(',',$sample)) {
			my $overlap;
			if ($num_samples == 1) {
				$overlap = 'no-overlap';
			} else {
				$overlap = 'overlap';
			}
			$results{'exp_concordance'}{$comparison}{$S}{$type}{$overlap} += $count;
			$results{'exp_concordance'}{$comparison}{$S}{'total'}{$overlap} += $count;
		}
		
		### add to sequence concordance
		next unless $num_samples == 2;
		my ($concordance);
		if ($num_genotypes == 1) {
			$concordance = 'concordant';
		} elsif ($num_genotypes == 2) {
			$concordance = 'discordant';
		} else {
			die "[STDERR]: more genotypes than expected\n";
		}
		
		$results{'seq_concordance'}{$comparison}{$type}{$zygosity}{$concordance} += $count;
		next if $zygosity eq 'het,homo-biallelic';
		### sum up by zygosity
		$results{'seq_concordance'}{$comparison}{$type}{'total'}{$concordance} += $count;
		### sum up by zygosity and type
		$results{'seq_concordance'}{$comparison}{'total'}{'total'}{$concordance} += $count; 
	}
}

####### print out results

foreach my $comparison (sort keys %{$results{'exp_concordance'}}) {
	foreach my $S (sort {$a<=>$b} keys %{$results{'exp_concordance'}{$comparison}}) {
		foreach my $type (sort keys %{$results{'exp_concordance'}{$comparison}{$S}}) {
			my @line = ('exp_concordance',$comparison,$relation{$comparison},$S,$other{$comparison}{$S},$type);
			foreach my $overlap ('overlap','no-overlap') {
				if (defined $results{'exp_concordance'}{$comparison}{$S}{$type}{$overlap}) {
					push(@line,$results{'exp_concordance'}{$comparison}{$S}{$type}{$overlap});
				} else {
					push(@line,0);
				}
			}
			print join("\t",@line), "\n";
		}
	}
}

foreach my $comparison (sort keys %{$results{'seq_concordance'}}) {
	foreach my $type (sort keys %{$results{'seq_concordance'}{$comparison}}) {
		foreach my $zygosity (sort keys %{$results{'seq_concordance'}{$comparison}{$type}}) {
			my @line = ('seq_concordance',$comparison,$relation{$comparison},$type,$zygosity);
			foreach my $concordance ('concordant','discordant') {
				if (defined $results{'seq_concordance'}{$comparison}{$type}{$zygosity}{$concordance}) {
					push(@line,$results{'seq_concordance'}{$comparison}{$type}{$zygosity}{$concordance});
				} else {
					push(@line,0);
				}
			}
			print join("\t",@line), "\n";
		}
	
	}
}
######################################################################################################

sub calc_type {
	my $genotypes = shift;
	
	### if any allele is greater than 1 in length, it's an insertion		
	my $type = 'snp';
	foreach my $genotype1 (split(',',$genotypes)) {
		foreach my $genotype (split(';',$genotype1)) {
			my @alleles = split('\|',$genotype);
			if (length($alleles[0]) > 1 || length($alleles[1]) > 1) {
				$type = 'insertion';
			}
		}
	}
	return $type;
}
