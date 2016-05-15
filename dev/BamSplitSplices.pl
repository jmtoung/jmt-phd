#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev', '/gpfs/fs121/h/toung/dev';
use File::Basename;

my $bam;

my $result = GetOptions( "bam=s" => \$bam );

my ($bam_name,$bam_dir,$bam_ext) = fileparse($bam,'\.bam');

$bam_ext eq '.bam' or die "[STDERR]: not a proper bam file: $bam\n";

open(COMMAND,"samtools view -h $bam |") or die "[STDERR]: can't open bamfile $bam with samtools\n";

my $sp_id = 1; ### id that groups together the splitted parts of an alignment

while(<COMMAND>) {
	if (/^@/) { print $_; next; } ### print the header

	chomp;
	my @split = split('\t');

	my $start = $split[3];
	my $cigar = $split[5];
	my $sequence = $split[9];
	my $qual = $split[10];
	
	if ($cigar =~ /N/) {
		my @cigar = split('(M|S|N|I|D|H)',$cigar);
		my @sequence = split('',$sequence);
		my @qual = split('',$qual);

		### these are all the new values. they get constantly updated, but first initialize them
		my $new_start = $start;
		my ($new_cigar,$new_sequence,$new_qual);
		
		my $start_shift = 0; ### $start_shift is how much you need to shift the start by
		my $last = 0; ### last says if you're on last iteration
		
		### if @cigar still defined, means you still have sequence to parse
		### if $new_sequence is still defined, means you still have sequence to print
		while(@cigar || $new_sequence) { 
			### check if @cigar is empty. if it is, that means we're on the last one...
			$last = 1 if !@cigar;

			my ($num,$type) = splice(@cigar,0,2);

			### if on last iteration or if type is 'N', we print a new record!
			### reason why we need $last variable is because we need to know if @cigar was empty PRIOR to splicing it. (Basically it takes care of the last iteration)
			if ($last || $type eq 'N') {
				my @new_split = @split; ### copy of alignment

				$new_split[3] = $new_start; ### set start equal to $new_start
				$new_split[5] = $new_cigar;
				$new_split[9] = $new_sequence;
				$new_split[10] = $new_qual;

				push(@new_split,'SP:i:'.$sp_id);
				push(@new_split,'ZR:Z:'.$start); ### ZR = original start
				push(@new_split,'ZC:Z:'.$cigar); ### ZC = original cigar
				push(@new_split,'ZS:Z:'.$sequence); ### ZS = original sequence
				push(@new_split,'ZQ:Z:'.$qual); ### ZQ = original quality
				print join("\t",@new_split), "\n"; ### print new record
				
				### shift the start by what you accumulated in $start_shift plus $num, which is the length of intron/gap
				$new_start += ($start_shift + $num) if $num; ### on the last iteration won't have $num defined
				$start_shift = 0; ### reset this since we already added it to $start
				$new_sequence = $new_qual = $new_cigar = ""; ### reset these values for next iteration

			} elsif ($type eq 'M' || $type eq 'S' || $type eq 'I' || $type eq 'D') {
				$start_shift += $num unless $type eq 'I' or $type eq 'S'; ### if sequence is clipped, don't need to shift start site
				$new_cigar .= $num . $type; ### append stuff to the next/new cigar 
				$new_sequence .= join("",splice(@sequence,0,$num)) unless $type eq 'D';
				$new_qual .= join("",splice(@qual,0,$num)) unless $type eq 'D';
			} else {
				die "[STDERR]: weird cigar string $split[5]\n";
			}
		}
		
		### check the @sequence and @qual are not defined
		!@sequence && !@qual or die "[STDERR] sequence and qual are still defined\n";
	} else {
		print $_, "\n";
	}
	
	$sp_id++;
}
close(COMMAND);
