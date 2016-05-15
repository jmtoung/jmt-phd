#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

my $gtf;
my $tag;
my $convertY2X;

my $options = GetOptions("gtf=s" => \$gtf, "tag=s" => \$tag, "convertY2X" => \$convertY2X);

######################################################################################
my ($gtf_name, $gtf_dir, $gtf_ext) = fileparse($gtf,'\.gtf');
### check that $FASTQ ends in '.fasta' and is an absolute file name
$gtf_ext eq '.gtf' or die "[STDERR]: '$gtf' does not end in '.gtf'\n";
substr($gtf_dir,0,1) eq '/' or die "[STDERR]: '$gtf' is not absolute file name\n";
### create output file
my $output = $gtf_dir . $gtf_name . '_' . $tag . $gtf_ext;
open(OUTPUT,'>'.$output) or die "[STDERR]: cannot open '$output': $!\n";
######################################################################################

open(GTF,$gtf) or die "[STDERR]: cannot open '$gtf': $!\n";
while(<GTF>) {
	chomp;
	my @split = split('\t');
	
	next if ($split[0] =~ /random|hap/);

	if ($split[0] eq 'Y' && $convertY2X) {
		my $region1_start = 1;
		my $region1_end = 2709520;
		my $region2_start = 57443438;
		my $region2_end = 57772954;
		my $region2_startX = 154584238;
		### if it's in the first region (1-2709520)
		if ($split[3] <= $region1_end && $split[4] <= $region1_end) {
			next;			
			#$split[0] = 'X'; 
		### if it's in the second region (57443438-57772954)
		} elsif ($split[3] <= $region2_end && $split[4] >= $region2_start) {
			next;			
			#$split[0] = 'X';
			#$split[3] = $split[3] - $region2_start + $region2_startX;
			#$split[4] = $split[4] - $region2_start + $region2_startX;
		### if it's only partially in the first region
		} elsif ($split[3] <= $region1_end && $split[4] > $region1_end) {
			print STDERR "[STDERR]: only_partially_in_first_region: $_\n";
		### if it's only partially in second region
		} elsif ($split[3] < $region2_start && $split[4] >= $region2_start || $split[3] <= $region2_end && $split[4] > $region2_end) {
			print STDERR "[STDERR]: only_partially_in_second_region: $_\n";
		}
	}	

	$split[0] =~ s/chr//;

	print OUTPUT join("\t",@split), "\n";
}
