#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use lib '/home/jmtoung/Lab/dev', '/gpfs/fs121/h/toung/oldhome/dev';

my $files;

my $options = GetOptions(
	"files=s" => \$files,
);


my @files = split(',',$files);

my $chrom;
foreach my $file (@files) {
	-e $file or die "[STDERR]: $file doesn't exist: $!\n";
	my ($file_name,$file_dir,$file_ext) = fileparse($file,'');
	my @file_name = split('\.',$file_name);
	my $file_chrom = $file_name[5];
	
	if (!defined $chrom) { $chrom = $file_chrom; }
	else { $chrom eq $file_chrom or die "[STDERR]: chromosomes differ: $chrom vs. $file_chrom\n"; }		
}

my $cat = "cat " . join(" ",@files) . " | uniq | sort +0 -1 +1n -2 +3n -4 +2 -3 |";
print STDERR $cat, "\n";
open(CAT,$cat) or die "[STDERR]: can't fork $cat: $!\n";
### %current holds current group of genotypes at site
### %results contains tallied results
### %samples contains all the samples in @files. to be used in output
my (%current,%results,%samples);
while(<CAT>) {
	chomp;
	
	my @split = split('\t');
	unless (is_equal(\@split,\%current)) {
		calculate_concordance(\%current);
		%current = ();
	} 
	
	add_to_current(\@split,\%current);
}

calculate_concordance(\%current);

### name the "comparison" by the sorted order of all samples we're comparing
my $comparison = join(",",sort {$a<=>$b} keys %samples);

foreach my $num_samples (sort {$a<=>$b} keys %results) {
	foreach my $samples (sort keys %{$results{$num_samples}}) {
		foreach my $num_genotypes (sort {$a<=>$b} keys %{$results{$num_samples}{$samples}}) {
			foreach my $genotypes (sort keys %{$results{$num_samples}{$samples}{$num_genotypes}}) {
				foreach my $zygosity (sort keys %{$results{$num_samples}{$samples}{$num_genotypes}{$genotypes}}) {
					my @results = ($comparison,$chrom,$num_samples,$samples,$num_genotypes,$genotypes,$zygosity,$results{$num_samples}{$samples}{$num_genotypes}{$genotypes}{$zygosity});
					print join("\t",@results), "\n";
				}
			
			}
		
		}
	
	}

}


######################################################################################

sub calculate_concordance {
	my $current = shift;
	
	### clean up redundant genotypes, e.g. removes A|A if AG|AG or AA|A exists
	foreach my $sample (keys %{$current->{'geno'}}) {
		foreach my $geno (sort keys %{$current->{'geno'}{$sample}}) {
			my ($a1,$a2) = split('\|',$geno);
			next if length($a1) == 1 && length($a2) == 1; ### next if both are single length alleles
			my $delete = substr($a1,0,1) . "|" . substr($a2,0,1); ### this is the thing to delete!, if AG|A, thing to delete is A|A
			delete $current->{'geno'}{$sample}{$delete};
		}
	}

	my (%samples,%genotypes,%zygosity);
	### calculate genotypes and zygosity for each sample
	foreach my $sample (keys %{$current->{'geno'}}) {
		### geno are all the genotype calls for this individual at this site.
		my @geno = sort keys %{$current->{'geno'}{$sample}};
			
		$samples{$sample}++;
		$genotypes{join(";",@geno)}++; ### amalgamate all genotypes and separate by ';'
		$zygosity{calc_zygosity(\@geno)}++; ### zygosity is either 'homo','het', or 'multi' (if there are more than 1 genotypes at this site)
	}
	
	my @samples = sort {$a<=>$b} keys %samples;
	my @genotypes = sort keys %genotypes;
	my @zygosity = sort keys %zygosity;

	### need to flag the special situation when zyogisty is het,homo.
	### need to distinguish between A|A| & A|G which could be due to low sampling/DAE
	### and A|A & C|T where clearly they are just discordant
	my $zygosity = join(",",@zygosity);
	if ($zygosity eq 'het,homo') {
		my %alleles;
		foreach my $genotype (@genotypes) {
			my ($a1,$a2) = split('\|',$genotype);
			$alleles{$a1}++;
			$alleles{$a2}++;
		}
		if (keys %alleles == 2) {
			$zygosity .= '-biallelic';
		}
	}
	
	$results{scalar(@samples)}{join(",",@samples)}{scalar(@genotypes)}{join(",",@genotypes)}{join(",",@zygosity)}++;
}

sub calc_zygosity {
	my $geno = shift;
	
	return 'multi' if @{$geno} > 1;
	
	my ($a1,$a2) = split('\|',$geno->[0]);
	
	return 'het' if $a1 ne $a2;
	return 'homo';
}

sub is_equal {
	my ($split,$current) = @_;
	
	return 1 if !defined $current->{'chrom'};
	
	return 0 if $current->{'chrom'} ne $split->[0];
	return 0 if $current->{'position'} != $split->[1];

	return 1;
}

sub add_to_current {
	my ($split,$current) = @_;
	
	$current->{'chrom'} = $split->[0];
	$current->{'position'} = $split->[1];
	$current->{'geno'}{$split->[3]}{$split->[2]}++;

	### add to the %samples hash to keep track of how many samples we have (to make the "comparison" tag)
	$samples{$split->[3]}++;	
}
