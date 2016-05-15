#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use lib '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use ComplementBase;

my $rdd_file;
my $blat_file;
my $rm_col;

my $options = GetOptions(
	"rdd_file=s" => \$rdd_file,
	"blat_file=s" => \$blat_file,
	"rm_col=s" => \$rm_col
);

print STDERR "rdd_file:\t$rdd_file\n";
print STDERR "blat_file:\t$blat_file\n";
print STDERR "rm_col:\t$rm_col\n";

my $num_sites_loaded;
my %rm_sites;
open(BLAT,$blat_file) or die "[STDERR]: can't open $blat_file: $!\n";
while(<BLAT>) {
	chomp;
	my @split = split('\t');
	
	foreach my $col (split(',',$rm_col)) {
		if ($split[$col] != 0) {
			### chrom, position, ref_base, rdd_base
			$num_sites_loaded++;
			$rm_sites{$split[0]}{$split[1]}{$split[3]}{$split[5]}++;
			goto NEXT;
		}
	}
	NEXT:
}
print STDERR "num_blat_sites_loaded:\t$num_sites_loaded\n";

my ($rdd_file_name,$rdd_file_dir,$rdd_file_ext) = fileparse($rdd_file,'\.(\w)+');
my $output = $rdd_file_dir . $rdd_file_name . "_rmblat" . $rdd_file_ext;
open(OUTPUT,">$output") or die "[STDERR]: can't open $output: $!\n";

open(RDD_FILE,$rdd_file) or die "[STDERR]: can't open $rdd_file: $!\n";
my $num_removed;
while(<RDD_FILE>) {
	chomp;
	my @split = split('\t');
	
	my $ref_base = "def";
	my $rdd_base = "def";
	if ($split[2] eq '-') { $ref_base = complement_base($split[3]); $rdd_base = complement_base($split[4]); }
	elsif ($split[2] eq '+' or $split[2] eq '.') { $ref_base = $split[3]; $rdd_base = $split[4]; }
	
	if (exists $rm_sites{$split[0]}{$split[1]}{$ref_base}{$rdd_base}) {
		$num_removed++;
	} else {
		print OUTPUT $_, "\n";
	}
}

print STDERR "num_removed\t$num_removed\n";

