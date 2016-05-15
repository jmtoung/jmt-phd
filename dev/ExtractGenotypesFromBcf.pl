#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use lib '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Database;

my $home;
my $bcf;
my $min_depth;
my $alnDB = "/ifs/h/toung/aln/alnDB";

my $options = GetOptions(
	"home=s" => \$home,
	"bcf=s" => \$bcf,
	"min_depth=s" => \$min_depth,
	"alnDB=s" => \$alnDB
);

defined $min_depth or die "[STDERR]: define min depth\n";
-e $bcf or die "[STDERR]: bcf $bcf doesn't exist: $!\n";

my ($bcf_name,$bcf_dir,$bcf_ext) = fileparse($bcf,'\.bcf');
my @bcf_name = split('\.',$bcf_name);

my $bam = $bcf_dir . join('.',@bcf_name[0..2]) . ".bam";
my $home_regex = $home; 
$home_regex =~ s/\//\\\//g;
$bam =~ s/$home_regex//g;
$bam =~ s/^\///;
my $alnID = Database->new($alnDB)->lookup(1,$bam,0);
defined $alnID or die "[STDERR]: can't find alnID for $bam in $alnDB\n";

my $output_dir = $bcf;
$output_dir =~ s/\.bcf/\.genotypes/;
unless (-d $output_dir) {
	mkdir($output_dir) or die "[STDERR]: can't make $output_dir directory\n";
}

my $below_depth = 0;
my $bcftools_view = "bcftools view $bcf |";
print $bcftools_view, "\n";
open(BCFTOOLS,$bcftools_view) or die "[STDERR]: can't open $bcftools_view: $!\n";
my %fh;
while(<BCFTOOLS>) {
	chomp;
	
	next if /^#/;
	
	my @split = split('\t');
	
	my $chrom = $split[0];
	my $pos = $split[1];
	my @alleles = ($split[3],split(',',$split[4]));
	my %info;
	foreach my $info (split(';',$split[7])) {
		my @info = split('=',$info);
		if (@info == 2) { $info{$info[0]} = $info[1]; }
		elsif (@info == 1) { $info{$info[0]}++; }
		elsif (@info >2) { die "[STDERR]: invalid format for info: $split[7]\n"; }
	}
	my %values;
	my @format = split(':',$split[8]);
	my @values = split(':',$split[9]);
	for (my $i=0; $i < @format; $i++) { $values{$format[$i]} = $values[$i]; }
	
	unless (defined $info{'DP4'} && get_depth($info{'DP4'}) >= $min_depth) {
		$below_depth++; next;
	}
	
	my $genotype; 
	if ($alleles[1] eq '.') {
		$genotype = $alleles[0] . "|" . $alleles[0];
	} else {
		$genotype = get_genotype(\@alleles,$values{'GT'});
	}

	unless (defined $fh{$chrom}) {
		my $output = $output_dir . "/" . $bcf_name . ".$chrom" . ".genotypes";
		open($fh{$chrom},">".$output) or die "[STDERR]: can't open $output: $!\n";
	}
	print {$fh{$chrom}} $chrom, "\t", $pos, "\t", $genotype, "\t", $alnID, "\n";
}

print "filtered-below_depth_$min_depth:\t$below_depth\t$bcf\n";

sub get_depth {
	my $DP4 = shift;
	
	my $depth = 0;
	foreach my $count (split(',',$DP4)) {
		$depth += $count;
	}
	return $depth;
}

sub get_genotype {
	my ($alleles,$GT) = @_;
	my @genotype;
	my @GT = split(/\/|\|/,$GT);
	foreach my $gt (@GT) {
		push(@genotype,$alleles->[$gt]);
	}
	
	return join("|",sort @genotype);
}
