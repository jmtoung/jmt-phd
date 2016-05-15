#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';

my $interval;

my $result = GetOptions(
	"interval=s" => \$interval
);

my ($start,$end);
if ($interval =~ /([0-9]+)-([0-9]+)/) {
	$start = $1;
	$end = $2;
} else {
	die "[STDERR]: 'interval needs to be num-num (e.g. 1-100)\n";
}

my %present;
while(<STDIN>) {
	chomp;
	$present{$_}++;
}

for(my $i = $start; $i <= $end; $i++) {
	print $i, "\n" unless exists $present{$i};
}

