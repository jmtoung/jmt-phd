#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

#my $file;
my $min_length = 35;
my $output;

my $options = GetOptions(
#	"file=s" => \$file,
	"min_length=i" => \$min_length,
	"output=s" => \$output
);

open(OUTPUT,">$output") or die "[STDERR]: can't open $output: $!\n";

#open(FILE,$file) or die "[STDERR]: cannot open '$file': $!\n";
#my ($file_name, $file_dir, $file_ext) = fileparse($file,'\.txt');

#### check that $FASTQ ends in '.fastq' and is an absolute file name
#$file_ext eq '.txt' or die "[STDERR]: '$file' does not end in '.txt'\n";
#substr($file_dir,0,1) eq '/' or die "[STDERR]: '$file' is not absolute file name\n";

#my $fastq = $file_dir . $file_name . ".fastq";
#open(FASTQ,'>'.$fastq) or die "[STDERR]: cannot open '$fastq': $!\n";
#select(FASTQ);

my $num_filter;
my $num_below_min_length;
while(my $line1 = <STDIN>) {
	defined (my $line2 = <STDIN>) or die "[STDERR]: line2 (sequence) undefined\n";
	defined (my $line3 = <STDIN>) or die "[STDERR]: line3 (read_name) undefined\n";
	defined (my $line4 = <STDIN>) or die "[STDERR]: line4 (qual scores) undefined\n";
	chomp ($line1, $line2, $line3, $line4);
	
	my @header = split('\s+',$line1);
	my @filter = split(':',$header[1]);
	
	if ($filter[1] eq 'Y') {
		$num_filter++;
		next;
	} elsif ($filter[1] eq 'N') {
	} else {
		die "[STDERR]: weird code\n";
	}
	
	my @line2 = split('',$line2);
	my @line4 = split('',$line4);

	### while @line4 is defined AND the last qual score is #, pop it off
	### @line4 might not be defined b/c the entire qual score might all be #'s
	### keep track of number popped off
	my $num_pop = 0;
	while(@line4 && $line4[-1] eq '#') {
		pop(@line4);
		pop(@line2);
		$num_pop++;
	}

	### if read(@line2) or quality score (@line4) are empty, next
	next unless (@line2 && @line4);
	
	### if read (@line2) length and quality score (@line4) length are not the same, exit
	if (scalar(@line2) != scalar(@line4)) {
		print STDERR "[STDERR]: 'sequence' and 'quality score' length unequal => $line1 | $line2 | $line3 | $line4\n";
		next;
	} 

	### if the length is not greater than min length, then skip
	if (scalar(@line2) < $min_length) {
		$num_below_min_length++;
		next;
	}

	$header[0] = $header[0] . "|$num_pop";
	my $header = join(' ',@header);
	print OUTPUT $header, "\n", join('',@line2), "\n", $line3, "\n", join('',@line4), "\n";
}

print "number filtered out ==> $num_filter\n";
print "number below min length ==> $num_below_min_length\n";
