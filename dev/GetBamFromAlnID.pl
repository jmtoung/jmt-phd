#!/usr/bin/perl -w

use strict;
use lib '/ifs/h/toung/dev', '/gpfs/fs121/h/toung/oldhome';
use Database;
use Getopt::Long;

my $home = "/gpfs/fs121/h/toung/oldhome";
my $alnDB = $home . "/aln/alnDB";

my $options = GetOptions(
	"home=s" => \$home,
	"alnDB=s" => \$alnDB
);

while(<STDIN>) {
	chomp;
	my @split = split('\t');
	print $split[0], "\t", Database->new($alnDB)->lookup(0,$split[0],1), "\n";
}
