#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $chrom_column;

my $options = GetOptions(
	"chrom_column=i" => \$chrom_column
);

$chrom_column =~ /[0-9]+/ or die "[STDERR]: chrom_column not valid\n";

while(<STDIN>) {
	chomp;
	my @split = split('\t');

	if ($split[$chrom_column] =~ s/^chr//) {
		if ($split[$chrom_column] =~ /^([0-9]{1,2}|X|Y|MT|M)$/) {
			$split[$chrom_column] = 'MT' if $split[$chrom_column] =~ /M/;
			print join("\t",@split), "\n";
		}
	}
}
