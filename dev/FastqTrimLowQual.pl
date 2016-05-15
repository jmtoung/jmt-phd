#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

my $fastq;
my $tag;
my $min_length;

my $result = GetOptions (
	"fastq=s" => \$fastq,
	"tag=s" => \$tag,
	"min_length=i" => \$min_length
);

my ($fastq_name, $fastq_dir, $fastq_ext) = fileparse($fastq,'\.fastq(\.gz)?');
my $fastq_fh;
if ($fastq_ext eq '.fastq') { $fastq_fh = $fastq; }
elsif ($fastq_ext eq '.fastq.gz') { $fastq_fh = "gunzip -c $fastq |"; }
else { die "[STDERR]: file doesn't end in '.fastq' or '.fastq.gz'"; }
open(FASTQ,$fastq_fh) or die "[STDERR]: cannot open '$fastq_fh': $!\n";

### check that $FASTQ ends in '.fastq' and is an absolute file name
#$fastq_ext eq '.fastq' or die "[STDERR]: '$fastq' does not end in '.fastq'\n";
substr($fastq_dir,0,1) eq '/' or die "[STDERR]: '$fastq' is not absolute file name\n";

my $fastq_new = $fastq_dir . $fastq_name . "_" . $tag . $fastq_ext;
$fastq_new =~ s/\.gz//;
open(FASTQ_NEW,'>'.$fastq_new) or die "[STDERR]: cannot open '$fastq_new': $!\n";
select(FASTQ_NEW);
my %num_trimmed;
my $num_removed;
while(my $line1 = <FASTQ>) {
	defined (my $line2 = <FASTQ>) or die "[STDERR]: line2 (sequence) undefined\n";
	defined (my $line3 = <FASTQ>) or die "[STDERR]: line3 (read_name) undefined\n";
	defined (my $line4 = <FASTQ>) or die "[STDERR]: line4 (qual scores) undefined\n";
	chomp ($line1, $line2, $line3, $line4);
	
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

	$num_trimmed{$num_pop}++;
	
	### if read(@line2) or quality score (@line4) are empty, next
	next unless (@line2 && @line4);
	
	### if read (@line2) length and quality score (@line4) length are not the same, exit
	if (scalar(@line2) != scalar(@line4)) {
		print STDERR "[STDERR]: 'sequence' and 'quality score' length unequal => $line1 | $line2 | $line3 | $line4\n";
		next;
	} 

	### if the length is not greater than min length, then skip
	if (scalar(@line2) < $min_length) {
		$num_removed++;
		next;
	}

	print $line1 . "|" . $num_pop, "\n", join('',@line2), "\n", $line3, "\n", join('',@line4), "\n";
}
close(FASTQ);

print STDOUT "num_removed:\t$num_removed\n";
print STDOUT "length\tcount\n";
foreach my $pop_length (sort {$a<=>$b} keys %num_trimmed) {
	print STDOUT $pop_length, "\t", $num_trimmed{$pop_length}, "\n";
}
