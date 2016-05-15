package CalculateBamStats;

require Exporter;
use strict;
use Getopt::Long;
use File::Basename;

our @ISA = qw(Exporter);
our @EXPORT = qw(calculate_bam_stats);

our ($TEMP_FILE);

sub calculate_bam_stats {
	my $BAM = shift;

	my ($bam_name, $bam_dir, $bam_ext) = fileparse($BAM,'\.bam');

	### check that $BAM ends in '.bam' and is an absolute file name
	$bam_ext eq '.bam' or die "[STDERR]: '$BAM' does not end in '.bam'\n";
	substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$BAM' is not absolute file name\n";

	### loop through bam file and gather some stats ########################
	my $PIPE = "samtools view $BAM | awk " . "'" . '$2!=4 {print $1}' . "' | sort | uniq -c |";
	open(PIPE,$PIPE) or die "[STDERR]: cannot fork: $!\n";
	my %hist;
	my $num_aligned;
	while(<PIPE>) {
		chomp;
		my @split = split('\s+');
		$hist{$split[1]}++;
		$num_aligned++;
	}

	my @hist;
	foreach my $times_aligned (sort {$a <=> $b} keys %hist) {
		push(@hist,$times_aligned.':'.$hist{$times_aligned});
	}
	my $hist = join(",",@hist);

	return $num_aligned,$hist;
}

END {
	### delete temporary files
	system("rm $TEMP_FILE") if ($TEMP_FILE && -e $TEMP_FILE);
}

1;
