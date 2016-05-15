#!/usr/bin/perl -w

use strict;
use lib '/gpfs/fs121/h/toung/dev';
use GetLineCount;
use Getopt::Long;
use File::Basename;

my $file;
my $lineSize;
my $outDir;

my $options = GetOptions(
	"file=s" => \$file,
	"lineSize=i" => \$lineSize,
	"outDir=s" => \$outDir
);

my ($fileName,$fileDir,$fileExt) = fileparse($file,'');

my $part = 1;
my $count = 0;
my %fh;

open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";
print STDERR "file:\t$file\n";

unless (defined $outDir) {
	$outDir = $fileDir;	
}

print STDERR "outDir:\t$outDir\n";

while(<FILE>) {
	if ($count == $lineSize) {
		$part++;
		$count = 0;
	}
	
	if (!defined $fh{$part}) {
		my $out = $outDir . "/" . $fileName . "." . formatInteger($part,7);
		open($fh{$part},">".$out) or die "[STDERR]: can't open $out: $!\n";
		if (defined $fh{$part - 1}) {
			close($fh{$part - 1});
		}
	}

	$count++;	
	print {$fh{$part}} $_;
}
close(FILE);

print STDERR "numPartsSplit:\t$part\t$file\n";
