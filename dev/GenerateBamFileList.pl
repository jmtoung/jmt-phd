#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Find;

my $directory;

my $options = GetOptions("directory=s" => \$directory, "output=s" => \$output);

my @FILES;

find(\&find_bam, $directory);

print join(" ",@FILES);

sub find_bam {
	next unless $_ =~ /\.bam$/;
	push(@FILES, $File::Find::name);
}
