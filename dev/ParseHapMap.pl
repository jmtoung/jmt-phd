#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use Cwd;

my $hapmap; # = "/home/jmtoung/Lab/rdd/hapmap/genotypes/2010-08_phaseII+III/forward/genotypes_chr10_CEU_r28_nr.b36_fwd.txt";
my $individuals; # = "NA12004,NA12750";

my $result = GetOptions(
	"hapmap=s" => \$hapmap,
	"individuals=s" => \$individuals
);

################################################################################
### This script extracts info from a HapMap file for individual(s) specified
################################################################################

################################################################################
my ($hapmap_name, $hapmap_dir, $ext) = fileparse($hapmap,'\.[a-zA-Z]+$');
### if $hapmap_dir is './', need to convert it to an absolute path
if ($hapmap_dir eq './') { $hapmap_dir = getcwd . '/'; }
my @individuals = split(',',$individuals);
################################################################################

################################################################################
### open hapmap file and get first line (header line)
open(HAPMAP,$hapmap) or die "couldn't fork: $!\n";
chomp (my $header_line = <HAPMAP>);
my @header_line = split('\s+',$header_line);
################################################################################

################################################################################
### get the indices of each individual and open output files for each
my %index; ### hash of the indices of each individual
my %output; ### hash of the output file handles for each individual
################################################################################

################################################################################
### find indices for individuals
foreach my $individual (@individuals) {

	### get the index for this individual
	my $index = 0;
	foreach (@header_line) {
		if ($_ eq $individual) { $index{$individual} = $index; }
		else { $index++; }
	}
	next unless (defined $index{$individual});

	########################################################################
	### in the hapmap_file directory, make a folder with the name of the hapmap file ($hapmap_name)
	unless(-d $hapmap_name) { mkdir $hapmap_name; }
	########################################################################

	### make a file in the hapmap_name dir for the individual's genotype/output
	chdir ($hapmap_dir);
	chdir ($hapmap_name);

	my $filename = $hapmap_name . "." . $individual . ".txt"; 
	open($output{$individual},'>'.$filename);
}
################################################################################

exit if scalar(keys %index) == 0;

while(<HAPMAP>) {
	chomp;
	my @line = split('\s+');

	foreach my $individual (@individuals) {
		my $chrom = $line[2];
		$chrom =~ s/chr//;
		my $pop_genotype = $
		### print chrom, position, strand, rs#, qc_code, genotype
		print {$output{$individual}} $line[2], "\t", $line[3], "\t", $line[4], "\t", $line[0], "\t", $line[10], "\t", $line[$index{$individual}], "\n";
	}
}
close(HAPMAP);
