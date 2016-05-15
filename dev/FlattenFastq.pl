#!/usr/bin/perl -w

use strict;

my $numReads = 0;

while(my $line1 = <STDIN>) {
	my $line2 = <STDIN>;
	my $line3 = <STDIN>;
	my $line4 = <STDIN>;
	
	chomp ($line1, $line2, $line3, $line4);
	
	print join("\t",$line2,$line1), "\n";
	
	$numReads++;
}

my $lineRemaining = <STDIN>;

!defined $lineRemaining or die "[STDERR]: still has lines\n";

print STDERR "read in $numReads fastq entries\n";
