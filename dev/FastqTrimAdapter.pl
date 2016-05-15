#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

my $fastq;
my $tag = "trimadapter";
my $adapter = "TCGTATGCCGTCTTCTGCTTG";
my $min_length = 35;
my $min_trim = 5;
my $no_mismatch_length = 10;
my $max_mismatch = 1;

my $options = GetOptions(
	"fastq=s" => \$fastq,
	"tag=s" => \$tag,
	"adapter=s" => \$adapter, 
	"min_length=i" => \$min_length,
	"min_trim=i" => \$min_trim,
	"no_mismatch_length=i" => \$no_mismatch_length,
	"max_mismatch=i" => \$max_mismatch
);

print STDOUT "fastq_file:\t$fastq\n";
print STDOUT "tag:\t$tag\n";
print STDOUT "adapter:\t$adapter\n";
print STDOUT "min_length:\t$min_length\n";
print STDOUT "min_trim:\t$min_trim\n";
print STDOUT "no_mismatch_length:\t$no_mismatch_length\n";
print STDOUT "max_mismatch:\t$max_mismatch\n";

my ($fastq_name, $fastq_dir, $fastq_ext) = fileparse($fastq,'\.fastq(\.gz)?');

my $fastq_command;
if ($fastq_ext eq '.fastq.gz') { $fastq_command = "gunzip -c $fastq |"; }
elsif ($fastq_ext eq '.fastq') { $fastq_command = $fastq; }
else { die "[STDERR]: what ext?\n"; }
open(FASTQ,$fastq_command) or die "[STDERR]: cannot open '$fastq_command': $!\n";

my $fastq_new = $fastq_dir . $fastq_name . "_" . $tag . $fastq_ext;
$fastq_new =~ s/\.gz//;
open(FASTQ_NEW,'>'.$fastq_new) or die "[STDERR]: cannot open '$fastq_new': $!\n";
select(FASTQ_NEW);

my %trim_stats;
while(my $line1 = <FASTQ>) {
	defined (my $line2 = <FASTQ>) or die "[STDERR]: line2 (sequence) undefined\n";
	defined (my $line3 = <FASTQ>) or die "[STDERR]: line3 (read_name) undefined\n";
	defined (my $line4 = <FASTQ>) or die "[STDERR]: line4 (qual scores) undefined\n";
	chomp ($line1, $line2, $line3, $line4);

	my $line2_preserve = $line2;

	### starting with the entire adapter, we get progressively smaller
	my $end_trimmed = 0;
	for (my $i = length($adapter); $i >= $min_trim; $i--) {
		### This is the piece of the adapter we're looking for in $read_seq
		my $adapter_piece = substr($adapter, 0, $i);
		### This is the end of the read_seq we believe contains the adapter
		my $read_piece = substr($line2, length($line2) - $i, $i);

		### if the two pieces match, then trim
		if (match($adapter_piece, $read_piece, $no_mismatch_length, $max_mismatch)) {
			### if length of what we're proposing to trim is greater than min_trim, then trim
			$end_trimmed = 1;
			
			### trim!
			substr($line2, length($line2) - $i, $i, '');
			substr($line4, length($line4) - $i, $i, '');
			
			### add to read name the length we trimmed
			$line1 .= "|" . length($adapter_piece);			
			last;
		}
	}
	
	### search for adapter in the entire read...but don't allow any mismatches
	my $mid_trimmed = 0;
	my $rindex = rindex($line2_preserve,$adapter); ### search in the full sequence
	unless($rindex == -1) {
		$mid_trimmed = 1;
		### if end is already trimmed, just make a note of it...
		unless ($end_trimmed) {
			my $trimmed_seq = substr($line2, $rindex, length($line2) - $rindex); ### sequence we're whacking off

			### whack off sequence
			$line2 = substr($line2,0,$rindex);
			$line4 = substr($line4,0,$rindex);
		
			### add to read name the length we trimmed
			$line1 .= "|" . length($trimmed_seq);
		}
	}


	
	### check to see that the read is greater than a certain length	
	print $line1, "\n", $line2, "\n", $line3, "\n", $line4, "\n" unless (length($line2) < $min_length);
	$trim_stats{$end_trimmed}{$mid_trimmed}++ unless (length($line2) < $min_length);
}

print STDOUT "end_trimmed\tmid_trimmed\tcount\n";
foreach my $end_trimmed (sort {$a<=>$b} keys %trim_stats) {
	foreach my $mid_trimmed (sort {$a<=>$b} keys %{$trim_stats{$end_trimmed}}) {
		print STDOUT $end_trimmed, "\t", $mid_trimmed, "\t", $trim_stats{$end_trimmed}{$mid_trimmed}, "\n";
	}
}

### this sub returns 1 or 0 depending on if it matched
sub match {
	my @adapter_piece = split('',shift);
	my @read_piece = split('',shift);
	my $no_mismatch_length = shift;
	my $max_mismatch = shift;
	
	### num mismatch we've found so far
	my $num_mismatch = 0;		
	foreach (my $j = 0; $j <= $#adapter_piece; $j++) {
		$num_mismatch++ if ($adapter_piece[$j] ne $read_piece[$j]);

		### if the read_piece is greater than the no_mismatch_length
		if(scalar(@adapter_piece) >= $no_mismatch_length) {
			### then, we allow $max_mismatch mismatches
			return 0 if ($num_mismatch > $max_mismatch);
		} else {
			return 0 if ($num_mismatch != 0);
		}
	}
	return 1;
}
