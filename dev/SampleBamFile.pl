#!/usr/bin/perl -w

use strict;
use File::Basename;
use Getopt::Long;

my $bam;
my $min_number = 10000000;
my $decay_rate = .5;

my $options = GetOptions(
	"bam=s" => \$bam,
	"min_number=i" => \$min_number,
	"decay_rate=s" => \$decay_rate
);

my ($bam_name,$bam_dir,$bam_ext) = fileparse($bam,'\.bam');

my $current_bam = $bam;
my $current_wc = get_wc($bam);

print "samplingDB", "\t", $current_bam, "\t", $current_wc, "\n";
my $count = 1;
do {
	
	($current_bam,$current_wc) = create_sampling($current_bam);
	print "samplingDB", "\t", $current_bam, "\t", $current_wc, "\n";
	$count++;
} until ($current_wc < $min_number);

sub create_sampling {
	my $current_bam = shift;

	my $new_bam = $bam_dir . $bam_name . ".split-$count" . $bam_ext;
	open(NEW_BAM,"| samtools view -bS - > $new_bam") or die "[STDERR]: can't fork $new_bam\n";
	my $command = "samtools view -h $current_bam |";
	print STDERR $command, "\n";
	open(BAM,$command) or die "[STDERR]: can't fork $command: $!\n";
	my $wc = 0;
	while(<BAM>) {
		if (/^@/) { print NEW_BAM $_; next; }
		
		my $rand = rand();
		if ($rand < $decay_rate) {
			$wc++;
			print NEW_BAM $_;
		}
	}
	close(BAM);

	return ($new_bam,$wc);
}

sub get_wc {
	my $bam = shift;
	
	my $wc = `samtools view $bam | wc -l`;
	chomp $wc;
	return $wc;	
}
