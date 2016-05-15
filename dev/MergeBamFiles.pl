#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $bam_files;

my $options = GetOptions(
	"bam_files=s" => \$bam_files
);

my @bam_files = split(',',$bam_files);

### check all bam_files exist ########################################
foreach my $bam_file (@bam_files) {
	-e $bam_file or die "[STDERR]: $bam_file doesn't exist: $!\n";
}
######################################################################

### take first bam file and output header ############################
my $first_bam_file = shift(@bam_files);
my $open_first_bam_file = "samtools view -h $first_bam_file |";
open(FIRST_BAM_FILE,$open_first_bam_file) or die "[STDERR]: can't fork\n";
my $num_lines = 0;
my $previous = undef;
while(<FIRST_BAM_FILE>) {
	if (/^@/) {
		print $_;
		next;
	}
	
	chomp;
	my @split = split('\t');
	
	unless (is_equal($split[0],$previous)) {
		$num_lines++;
		$previous = $split[0];
	}
	
	print $_, "\n";
	
}
print STDERR join("\t",'num_lines',$first_bam_file,$num_lines), "\n";
######################################################################

### open all the others... ###########################################
foreach my $bam_file (@bam_files) {
	my $open_bam_file = "samtools view $bam_file |";
	open(BAM_FILE,$open_bam_file) or die "[STDERR]: can't fork\n";
	$num_lines = 0;
	$previous = undef;
	while(<BAM_FILE>) {
		chomp;
		my @split = split('\t');
		
		unless (is_equal($split[0],$previous)) {
			$num_lines++;
			$previous = $split[0];
		}
		print $_, "\n";
		
	}
	print STDERR join("\t",'num_lines',$bam_file,$num_lines), "\n";
}
######################################################################

sub is_equal {
	my ($readname,$previous) = @_;
	
	return 0 if !defined $previous;
	
	return 0 if $previous ne $readname;
	
	return 1;
}
