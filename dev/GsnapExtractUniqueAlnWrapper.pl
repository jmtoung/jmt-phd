#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use File::Basename;

my $bam;
my $min_mapqual;

my $options = GetOptions(
	"bam=s" => \$bam,
	"min_mapqual=i" => \$min_mapqual
);

################################################################################
### This script produces unique alignments from Gsnap bam file
################################################################################

### LOAD $bam AND MAKE OUTPUT ##################################################
my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
$bam_ext eq '.bam' or die "[STDERR]: '$bam' does not end in '.bam'\n";
substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$bam' is not absolute file name\n";
################################################################################

my $samtools = "samtools view -h $bam |";
open(SAMTOOLS,$samtools) or die "[STDERR]: can't fork $samtools: $!n";

my $current_split0;
my @current;

while(<SAMTOOLS>) {
	if (/^@/) {
		print $_;
		next;
	}
	chomp;
	my @split = split('\t');
	if (!defined $current_split0 || $split[0] ne $current_split0) {
		process_current(\@current);
		@current = ();
		$current_split0 = $split[0];
	}
	push(@current,\@split);
}

process_current(\@current);

sub process_current {
	my $current = shift;

	### check that there's only one alignment
	if (scalar(@{$current}) == 1) {
		### check that flag is 0 or 16
		my $flag = $current->[0]->[1];
		return unless $flag == 0 or $flag == 16;
	
		### check that mapqual is greater than minqual
		my $mapqual = $current->[0]->[4];
		return unless $mapqual >= $min_mapqual;
		
		print join("\t",@{$current->[0]}), "\n";
	}	
}
