#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use File::Basename;

my $bam;

$|++;
my $result = GetOptions(
	"bam=s" => \$bam
);

my ($bam_name,$bam_dir,$bam_ext) = fileparse($bam,'\.bam');

my $command = "samtools view -h $bam |";
print STDERR $command, "\n";
open(COMMAND,$command) or die "[STDERR]: can't fork $command: $!\n";

while(<COMMAND>) {
	if (/^@/) { print $_; next; }

	chomp;
	my @split = split('\t');

	my $start = $split[3];
	my $cigar = $split[5];
	my @sequence = split('',$split[9]);
	my @qual = split('',$split[10]);
		
	if ($cigar =~ /N/) {
		
		my @cigar = split('(M|S|N|I|D)',$cigar);

		my ($accum_bases,$last);
		my ($cigar,$sequence,$qual); ### new values
		while(@cigar || $sequence) {
			### check if @cigar is empty. if it is, that means we're on the last one...
			$last = 1 if (scalar(@cigar) == 0);

			my ($num,$type) = splice(@cigar,0,2);

			if ($last || $type eq 'N') {
				my @new_split = @split;
				$new_split[3] = $start;
				$new_split[5] = $cigar;
				$new_split[9] = $sequence;
				$new_split[10] = $qual;
				print join("\t",@new_split), "\n";

				$start += $accum_bases + $num if $num; ### on the last iteration won't have $num defined
				$sequence = $qual = $cigar = "";			
			} elsif ($type eq 'M' || $type eq 'S' || $type eq 'I' || $type eq 'D') {
				$cigar .= $num . $type;
				$sequence .= join("",splice(@sequence,0,$num)) unless $type eq 'D';
				$qual .= join("",splice(@qual,0,$num)) unless $type eq 'D';
				$accum_bases += $num unless $type eq 'I';
			} else {
				die "[STDERR]: weird cigar string $split[5]\n";
			}
		}
	} else {
		print $_, "\n";
	}	
}
close(COMMAND);
