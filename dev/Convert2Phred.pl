#!/usr/bin/perl -w

use strict;
use lib '/gpfs/fs121/h/toung/dev';
use ConvertAscii2Phred;

while(<STDIN>) {
	chomp;
	
	my @split = split('\t');
	
	my $qual = convert_quality($split[5]);
	
	print join("\t",@split,join(",",@{$qual})), "\n";

}

