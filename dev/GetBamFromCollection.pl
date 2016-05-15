#!/usr/bin/perl -w

use strict;
use lib '/ifs/h/toung/dev', '/gpfs/fs121/h/toung/oldhome';
use Database;
use Getopt::Long;

my $home = "/gpfs/fs121/h/toung/oldhome";
my $collection;
my $alnDB = "$home/aln/alnDB";

my $options = GetOptions(
	"home=s" => \$home,
	"collection=s" => \$collection,
	"alnDB=s" => \$alnDB
);

$collection = "$home/oldhome/database/collectionsDB/" . $collection;

open(COLLECTION,$collection) or die "[STDERR]: can't open $collection:$! \n";

while(<COLLECTION>) {
	chomp;
	my @split = split('\t');
	print $split[1], "\t", Database->new($alnDB)->lookup(0,$split[1],1), "\n";
}
