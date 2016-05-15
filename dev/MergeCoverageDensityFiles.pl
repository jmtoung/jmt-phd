#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Find;

my $cov_dirs;

my $options = GetOptions(
	"cov_dirs=s" => \$cov_dirs
);

### HASH THAT STORES COVERAGE INFO
my %COV;

### LOOP THROUGH EACH DIRECTORY
foreach my $cov_dir (split(',',$cov_dirs)) {

	opendir(COV_DIR,$cov_dir) or die "[STDERR]: can't open $cov_dir: $!\n";
	my @files = readdir(COV_DIR);
	
	### FOR EACH FILE THAT ENDS IN '.cov'
	foreach my $file (@files) {
		next unless $file =~ /\.cov$/;
		my $full_file = $cov_dir . ((substr($cov_dir,-1,1) ne '/') && "/") . $file;
		open(FILE,$full_file) or die "[STDERR]: can't open $full_file: $!\n";
		while(<FILE>) {
			chomp;
			my ($coverage, $count, $aln_id) = split('\t');
			$COV{$aln_id}{$coverage} += $count;
		}
	}
}
###

### PRINT DATA
foreach my $aln_id (sort {$a <=> $b} keys %COV) {
	foreach my $coverage (sort {$a <=> $b} keys %{$COV{$aln_id}}) {
		print $coverage, ",", $COV{$aln_id}{$coverage}, ",", $aln_id, "\n";
	}
}
