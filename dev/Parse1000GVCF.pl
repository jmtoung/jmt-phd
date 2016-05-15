#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use Cwd;

my $vcf_file = "/home/jmtoung/Lab/rdd/1000_genomes/pilot_data/release/2010_07/exon/snps/test.txt";
my $identifier = "snps_PASS";
my $type = "SNP";
my $individuals = "06985,06986";
my $filter = 1;

my $result = GetOptions(
	"vcf_file=s" => \$vcf_file, 
	"identifier=s" => \$identifier, 
	"type=s" => \$type, 
	"individuals=s" => \$individuals,
	"filter=i" => \$filter	
);

################################################################################
### This script takes a VCF file and outputs a list of genotypes from VCF for individual(s) specified. 
################################################################################

################################################################################
my ($vcf_name, $vcf_dir, $vcf_ext) = fileparse($vcf_file,'\.[a-zA-Z]+$');
if ($vcf_dir eq './') { $vcf_dir = getcwd . '/'; }
my @individuals = split(',',$individuals);
################################################################################

################################################################################
### open vcf file and get first line (header line)
open(VCF_FILE, $vcf_file) or die "couldn't fork: $!\n";
################################################################################

################################################################################
### inside directory with the file, make folder with vcf_file name and the identifier name
chdir($vcf_dir);
my $new_dir = $vcf_name . "_" . $identifier;
unless (-d $new_dir) { mkdir $new_dir; }
chdir($new_dir);
################################################################################

################################################################################
my %index; ### hash of indices for each individual
my %output; ### hash of output filehandles for each individual
################################################################################
my @sample_line;
################################################################################
while(<VCF_FILE>) {
	chomp;
	my @line = split('\t');
	
	### these are info lines
	if (/^#/) {
		### this line contains the sample IDs
		if (/^#CHROM/) {
			@sample_line = @line;
			
			foreach my $individual (@individuals) {
				my $col = get_column(\@sample_line,("NA".$individual));
				$index{$individual} = $col if $col;	
			}
		}
		next;
	}
	
	### note only individuals in %index were found/exist in this file
	foreach my $individual (keys %index) {

		unless (exists $output{$individual}) {
			chdir($vcf_dir);
			chdir($new_dir);
			my $file = $vcf_name . "." . $identifier . ".NA" . $individual . ".vcf";
			open($output{$individual},'>'.$file);
		}
	
		### print to output file
		select($output{$individual});
		### get_snp(ref_base,alt_bases,format,raw_genotype)
		if ($type eq 'SNP') { 
			### require the length of reference allele to be 1
			next unless length($line[3]) == 1;
			### require the length of each alternative allele to be 1
			foreach my $alt (split(',',$line[4])) { goto NEXT if length($alt) != 1; }
			### if there is a filter...
			if ($filter) { next unless $line[6] eq 'PASS'; }

			### print stuff
			my $id = $line[2];
			$id = "thouG" if $id eq '.';
			###  chrom, position, genotype
			print $line[0], "\t", $line[1], "\t", $id, "\t", "SNP", "\t", get_snp($line[3],$line[4],$line[8],$line[$index{$individual}]), "\n";
		}
		### get_indel()
		if ($type eq 'indel') {
			my @format = split(':',$line[8]);
			my @genotype_raw = split(':',$line[$index{$individual}]);
			foreach (my $i = 0; $i <= $#format; $i++) {
				next unless $format[$i] eq 'GT';
				if ($genotype_raw[$i] =~ /1/) {
					my @ref = split('',$line[3]);
					my @var = split('',$line[4]);
					### this is an insertion
					if ($#ref < $#var) {
						print $line[0], "\t", $line[1], "\t", "INS", "\t", $line[4], "\n";
					} elsif ($#ref > $#var) {
						for (my $i = ($#var + 1); $i <= $#ref; $i++) {
							print $line[0], "\t", ($line[1] + $i), "\t", "DEL", "\t", "D", "\n";
						}
					} else {
						for (my $i = 0; $i <= $#var; $i++) {
							if ($ref[$i] ne $var[$i]) {
								print $line[0], "\t", ($line[1] + $i), "\t", "SUB", "\t", $var[$i], "\n";
							}
						}
					}
				}
			}
		}
		NEXT:
	}

}
################################################################################

sub get_column {
	my $line = shift;
	my $sample = shift;
	
	my $col = 0;
	
	foreach my $val (@{$line}) {
		return $col if $val eq $sample;
		$col++
	}
	return undef;
}

sub get_snp {
	my $ref = shift;
	my @alt = split(',',shift(@_)); ### there might be more than one alternative allele
	my @format = split(':',shift(@_));
	my @genotype_raw = split(':',shift(@_));
	
	### return this
	my $genotype;
	
	### find the 'GT' field
	foreach (my $i = 0; $i <= $#format; $i++) {
		next unless $format[$i] eq 'GT';
		my @GT_field = split('',$genotype_raw[$i]);
		while(@GT_field) {
			my $genotype_raw = shift(@GT_field);
			if ($genotype_raw =~ /^0$/) {
				$genotype .= $ref;
			} elsif ($genotype_raw =~ /^([0-9]+)$/) {
				$genotype .= $alt[$1-1];
			} elsif ($genotype_raw =~ /\./) {
				$genotype .= '.';
			}
		}
		return $genotype;
	}
}
