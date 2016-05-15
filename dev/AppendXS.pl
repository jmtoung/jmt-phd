#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use File::Basename;

my $bam;
my $ss_map;
my $region;

my $result = GetOptions(
	"bam=s" => \$bam,
	"ss_map=s" => \$ss_map,
	"region=s" => \$region
);

### LOAD BAM ###################################################################
my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
$bam_ext eq '.bam' or die "[STDERR]: '$bam' does not end in '.bam'\n";
substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$bam' is not absolute file name\n";
################################################################################

### MAKE OUTPUT DIRECTORY & FILES ##############################################
#my $output = $bam . "_xs"; 
################################################################################

################################################################################
my %SS_MAP;
open(SS_MAP,$ss_map) or die "[STDERR]: can't open $ss_map: $!\n";
while(<SS_MAP>) {
	chomp;
	next if (/^###/);
	my @split = split('\t');
	my @split2 = split(';',$split[2]);
	foreach my $split2 (@split2) {
		my @split3 = split(':',$split2);
		push(@{$SS_MAP{$split[0]}},$split3[0]);
	}
}

my $samtools_view = "samtools view -h $bam |";
### if we only want a region of the bam file
$samtools_view = "samtools view -h $bam $region |" if $region;
open(SAMTOOLS,$samtools_view) or die "[STDERR]: can't fork $samtools_view: $!\n";
while(<SAMTOOLS>) {
	chomp;
	if (/^@/) {
		print $_, "\n";
	} else {
		my @split = split('\t');
		if ($split[5] =~ /N/) {
			if ($_ =~ /SS:Z:(ss[0-9]+)/) {
				my @strands = $SS_MAP{$1};
				if (scalar(@strands) == 1) {
					if ($strands[0] eq '+') { 
						push(@split,"XS:A:+");
					} else {
						push(@split,"XS:A:-");
					}
					print join("\t",@split), "\n";	
				} else {
					print "[STDERR]: more than 1 
strand\n";
				}
				
			} else {
				print "[STDERR]: no ss # for $_\n";
			}
		} else {
			print $_, "\n";
		}
	}
}
