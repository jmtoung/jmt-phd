#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $interval;
my $count;
my $fileName = "sge_RunGsnap.pl";

my $options = GetOptions(
	"interval=s" => \$interval,
	"count=i" => \$count,
	"fileName=s" => \$fileName
);

print "*****REPORT******\n\n";

my ($start,$end) = split('-',$interval);

my %found;
foreach my $c ($count - 1,$count,$count + 1) {

	my $command = "grep '$c' $fileName.* |";

	open(COMMAND,$command) or die "[STDERR]: can't open $command: $!\n";

	while(<COMMAND>) {
		chomp;

		my @split = split(':');
	
		my @split2 = split('\.',$split[0]);
		
		$found{$c}{$split2[3]}++;
	} 
	close (COMMAND);
}

my $total_finished = 0;
my $total_lines_finished = 0;
foreach my $c (sort {$a<=>$b} keys %found) {
	my @finished = keys %{$found{$c}};
	my $numFinished = scalar(@finished);
	print $numFinished, " jobs finished with $c lines\n";
	$total_lines_finished += $numFinished * $c;
	$total_finished += $numFinished;
}

my $sep = "----------------------------------------";
print "$sep\n";

print $total_finished, " jobs finished in total\n";
print $total_lines_finished, " lines finished in total\n";

my $prev = -100;
for (my $i = $start; $i <= $end; $i++) {
	
	my $found = 0;
	
	foreach my $c (keys %found) {
		if (defined $found{$c}{$i}) {
			$found = 1;
		}
	}
	
	if (!$found) {
		unless ($i == $prev + 1) {
			print "$sep\n";
		}
		print $i, "\tnotFinished\n" if !$found;	
		$prev = $i;
	}
}
