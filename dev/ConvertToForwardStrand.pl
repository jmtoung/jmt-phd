#!/usr/bin/perl -w

use strict;
use lib '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use ComplementBase;

while(<STDIN>) {
	chomp;
	my @split = split('\t');

	if ($split[2] eq '+') { $split[2] = '+,-'; }
	if ($split[2] eq '-') {
		$split[2] = '+,-';
		$split[3] = complement_base($split[3]);
		$split[5] = complement_base($split[5]);
	}

	print join("\t",@split), "\n";
}
