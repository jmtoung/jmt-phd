#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use File::Basename;

my $bam;
my $max_splice_length = 2000000;

$|++;
my $result = GetOptions(
	"bam=s" => \$bam,
	"max_splice_length=i" => \$max_splice_length
);

my ($bam_name,$bam_dir,$bam_ext) = fileparse($bam,'\.bam');

my $command = "samtools view -h $bam |";
print STDERR $command, "\n";
open(COMMAND,$command) or die "[STDERR]: can't fork $command: $!\n";

my $num_spliced;
my $num_removed;
my %splice_sizes;
while(<COMMAND>) {
	if (/^@/) { print $_; next; }
	chomp;

	my @split = split('\t');
	my @cigar = split('(M|S|N|I|D)',$split[5]);
	my $total_splice_size = 0;
	while(@cigar) {
		my $num = shift(@cigar);
		my $type = shift(@cigar);
		if ($type eq 'N') { $total_splice_size += $num; }
	}

	if ($total_splice_size) {
		$num_spliced++;
		$splice_sizes{$total_splice_size}++;
		if ($total_splice_size > $max_splice_length) {
			$num_removed++;
			print STDERR "rm_aln", "\t", $_, "\n";
		} else {
			print $_, "\n";
		}
	} else {
		print $_, "\n";
	}
}
close(COMMAND);

print STDERR "num_spliced", "\t", $num_spliced, "\n";
print STDERR "num_removed", "\t", $num_removed, "\n";

foreach my $splice_size (sort {$a<=>$b} keys %splice_sizes) {
	print STDERR "splice_size_count", "\t", $splice_size, "\t", $splice_sizes{$splice_size}, "\n";
}

my %count_thresholds;
for (my $threshold = 10000; $threshold <= 500000; $threshold += 10000) {
	foreach my $splice_size (sort {$a<=>$b} keys %splice_sizes) {
		if ($splice_size <= $threshold) {
			$count_thresholds{$threshold} += $splice_sizes{$splice_size};
		}
	}
}

foreach my $threshold (sort {$a<=>$b} keys %count_thresholds) {
	print STDERR "num_lt_threshold", "\t", $threshold, "\t", $count_thresholds{$threshold}, "\n";
}
