#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $file;
my $minRddLevel; # = 5; 
my $minRddCount; # = 1; 
my $minTotalCount; # = 10;

my $options = GetOptions(
	"file=s" => \$file,
	"minRddLevel=s" => \$minRddLevel,
	"minRddCount=s" => \$minRddCount,
	"minTotalCount=s" => \$minTotalCount,
);

print STDERR "file:\t$file\n";
(defined $minRddLevel && print STDERR "minRddLevel:\t$minRddLevel\n") or die "[STDERR]: minRddLevel not defined\n";
(defined $minRddCount && print STDERR "minRddCount:\t$minRddCount\n") or die "[STDERR]: minRddCount not defined\n";
(defined $minTotalCount && print STDERR "minTotalCount:\t$minTotalCount\n") or die "[STDERR]: minTotalCount not defined\n";

my $FH = $file;
$FH = "zcat $file |" if $file =~ /\.gz$/;
open(FILE,$FH) or die "[STDERR]: can't open $file: $!\n";

while(<FILE>) {
	chomp;
	
	my @split = split('\t');
	
	### determine whether this site should be called RDD or not
	
	my $refBase = $split[4];
	my $totalCount = $split[6];

	### filter on total count
	next unless $totalCount >= $minTotalCount;
	
	my @bases = ('A','C','G','T');
	
	for (my $i = 0; $i < @bases; $i++) {
		next if $bases[$i] eq $refBase; ### skip if base is refBase b/c can't be rdd

		my $baseCount = $split[7 + $i];

		### filter on rddCount
		next unless $baseCount >= $minRddCount;
		
		### filter on rddLevel
		my $rddLevel = $baseCount/$totalCount*100;
		
		next unless $rddLevel >= $minRddLevel;
		
		my @line = (@split[0..4],$bases[$i],$totalCount,$baseCount,$.);
		print join("\t",@line), "\n";
	}
	
}

close(FILE);
