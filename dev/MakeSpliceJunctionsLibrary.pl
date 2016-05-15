#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use List::Util qw(max min);
use lib qw(
/ifs/h/toung/usr/local/lib/perl5
/ifs/h/toung/usr/local/lib/perl5/site_perl
);
use Bio::SeqIO;
use Bio::Seq;

my $output; ### name of output file
my $gtf; # = "/home/jmtoung/Lab/database/gencode.v3c.annotation.NCBI36.gtf,/home/jmtoung/Lab/database/refseqgenes_refgene_03072011.gtf"; 
my $index; # = "/home/jmtoung/Lab/database/hg18/hg18_all.fa";
my $read_length; # = 50;

my $options = GetOptions("output=s" => \$output, "gtf=s" => \$gtf, "index=s" => \$index, "read_length=s" => \$read_length);

######################################################################################
my ($output_name, $output_dir, $output_ext) = fileparse($output,'\.[^.]+');
### check that $FASTQ ends in '.fasta' and is an absolute file name
($output_ext eq '.fasta' || $output_ext eq '.fa') or die "[STDERR]: '$output' does not end in '.fasta' or '.fa'\n";
substr($output_dir,0,1) eq '/' or die "[STDERR]: '$output' is not absolute file name\n";
### create map file
my $output_map = $output_dir . $output_name . '.map';
open(MAP,'>'.$output_map) or die "[STDERR]: cannot open '$output_map': $!\n";
print MAP "### output: ", $output, "\n";
print MAP "### gtf: ", $gtf, "\n";
print MAP "### index: ", $index, "\n";
print MAP "### read_length: ", $read_length, "\n";
######################################################################################

######################################################################################
### open index fai file to read in chromosomes
my $index_fai = $index . ".fai";
open(INDEX_FAI,$index_fai) or die "[STDERR]: '$index_fai' does not exist\n";
my %chrom;
while(<INDEX_FAI>) {
	chomp;
	my @split = split('\t');
	$chrom{$split[0]}++;
}
######################################################################################

######################################################################################
### fasta file
my $fasta = Bio::SeqIO->new(-file => ">".$output, -format => 'fasta');
######################################################################################

######################################################################################
### hash that keeps track of all splice sites
my %splice_sites;
######################################################################################

######################################################################################
### put records into $gtf file
foreach my $GTF (split(',',$gtf)) {

	### hash that keeps the chrom, start, and end for all transcripts of this gtf file
	my %GTF;
	
	open(GTF,$GTF) or die "[STDERR]: cannot open '$GTF': $!\n";
	
	### loop through GTF and store chrom, starts and ends into hash
	while(<GTF>) {

		### skip header lines
		next if (/^##/);

		chomp;

		### split the line
		my @split = split('\t');
		
		### skip unless it's a chromosome in the index
		next unless (exists $chrom{$split[0]});

		### skip unless it's an exon
		next unless ($split[2] eq 'exon');
		
		### we are going to split column 9, the attributes column and store the keys and values in a hash
		my %split2;
		foreach (split(';',pop(@split))) {
			if (/\s*(.*)\s"(.*)"/) {
				$split2{$1} = $2;
			} elsif (/level\s([0-9])+/) {
				$split2{"level"} = $1;
			}
		}
		
		### store variables into hash
		my $index = $split2{"transcript_id"};
		my $chrom = $split[0];
		my $start = $split[3];
		my $end = $split[4];
		### store the number of times this chromosome is seen (to check later and make sure its equivalent to the number of exons).
		$GTF{$index}{$chrom}{"count"}++;
		### store the start and end
		$GTF{$index}{$chrom}{"start-end"}{$start}{$end}++;
	}
	
	### now that we have all the data from the gtf file stored in a hash, we can add to the splice sites library
	foreach my $index (keys %GTF) {
		
		### for each chromosome this transcript is on.
		foreach my $chrom (keys %{$GTF{$index}}) {
		
			### get all the exon starts and ends into arrays, easier to manipulate
			my (@starts, @ends);
			foreach my $start (sort {$a <=> $b} keys %{$GTF{$index}{$chrom}{"start-end"}}) {
				push(@starts,$start);
				foreach my $end (keys %{$GTF{$index}{$chrom}{"start-end"}{$start}}) {
					push(@ends,$end);
				}
			}

			### check that the size of the starts and ends array is the same as the number of times this chromosome was "seen" in the gtf file for this transcript
			if (scalar(@starts) != $GTF{$index}{$chrom}{"count"} || scalar(@ends) != $GTF{$index}{$chrom}{"count"}) {
				print STDERR "[STDERR]: error in number of exons: $index\n";
				next;
			} 

			### this hash keeps all the starts and ends of the pieces you're going to yoink out FOR THIS TRANSCRIPT!!!
			my %junction_boundaries;
		
			### loop through and add to splice junctions library
			### $i is the exon number you're on. 
			for (my $i = 0; $i < $#starts; $i++) {		

				### this is for the part of the read to the LEFT of the junction site
				### this index $j allows you to go to other exons relative to $i if need be.
				my $j = $i;
				### keep track of number of bases added to the left.
				my $bases_added = 0;
				while ($j >= 0 && $bases_added != $read_length - 1) {
					### the start of this junction read will be either the end of the current exon minus one less than the read length OR the start of this exon, whichever is greater
					my $junction_start = max($starts[$j], $ends[$j] - (($read_length - 1) - $bases_added) + 1);
					### add to the hash that keeep track of starts and ends for THIS transcript.
					### the start is what we just calculated. the end will always be the end of the exon you're on.
					$junction_boundaries{$i}{"left"}{$junction_start}{$ends[$j]}++;
					### update the number of bases added
					$bases_added += $ends[$j] - $junction_start + 1;
					### decrease the exon you're on. if the number of bases added is equal to one less the read length, then the loop will not execute again. only if the number of bases_added
					### is not equal to one less the read length will we go to the exon below. 
					### also, if $j is less than 0, loop also will not execute!
					$j--;
				}

				### this is for the part of the read to the RIGHT of the junction site
				### set index $j back to the exon we are on
				$j = $i;
				### keep track of number of bases added to the right
				$bases_added = 0;
				while ($j < $#starts && $bases_added != $read_length - 1) {
					### end of this junction read will either be the end of the next exon OR the start of the next exon plus whatever's remaining 
					my $junction_end = min($ends[$j + 1], $starts[$j + 1] + (($read_length - 1) - $bases_added) - 1);
					### add to the hash that keep track of starts and ends for THIS transcript.
					### the end is what we just calculated. the start will always be the start of the NEXT exon you're on.
					$junction_boundaries{$i}{"right"}{$starts[$j + 1]}{$junction_end}++;
					### update the number of bases added
					$bases_added += $junction_end - $starts[$j + 1] + 1;
					### increase the exon you're on. if the number of bases added is equal to one less than the read length, then the loop will not execute again.
					$j++;
				}
			}
		
			### add the entire series of starts and ends for this transcript to to the entire hash library
			### $i refers to the exon you're on (the $i variable from above).
			foreach my $i (sort {$a <=> $b} keys %junction_boundaries) {
				### these are the starts and ends for all sites to the left of the junction site
				my @left_boundaries;
				### these are the starts and ends for all sites to the right of the junction site
				my @right_boundaries;
				### for all the starts and ends to the left, add to the array
				foreach my $start (sort {$a <=> $b} keys %{$junction_boundaries{$i}{"left"}}) {
					foreach my $end (keys %{$junction_boundaries{$i}{"left"}{$start}}) {
						push(@left_boundaries,$start.",".$end);
					}
				}
				### for all starts and ends to the right, add to the array
				foreach my $start (sort {$a <=> $b} keys %{$junction_boundaries{$i}{"right"}}) {
					foreach my $end (keys %{$junction_boundaries{$i}{"right"}{$start}}) {
						push(@right_boundaries,$start.",".$end);
					}
				}
				
				### the key of the splice sites looks like this: chr1,chr4|2,6;7,29|40,67;80,90. The value is a hash with key chr1,chr4|2,6;7,29;40,67;80,90 and value being the transcript id.
				### These two strings are very similar. The difference is that the first one does not distinguish between sequences to the left and right of the splice site, whereas the second one does,
				### as signified by a |. THe first one is used to compare uniquess, because essentially these are the same sequences. They only differ in the location of the splice site.
				### the latter information is needed for downstream analyses.
				my $key = $chrom . "|" . join(';',@left_boundaries) . ";" . join(';',@right_boundaries);
				my $value = $chrom . "|" . join(';',@left_boundaries) . "|" . join(';',@right_boundaries);
				$splice_sites{$key}{$value}{$index}++;
			}
		}
	}
	close(GTF);
}

### now we loop through the splice sites hash and pull the sequences and print to output files
### this is the index that will be used to map sequences to transcript id, etc.
my $ss_index = 1;
### the $ss_location (or key of %splice_sites) is the $key from above, look up a couple of lines, this tells you the start and end of all parts to yoink out.
foreach my $ss_location (keys %splice_sites) {
	my @ss_location = split('\|',$ss_location);
	my $chrom = shift(@ss_location);
	
	### the splice_junction read, the actual read to be written into the fasta file!!!
	my $ss_read;
	### the starts and ends of each "chunk" or part of the read we are going to pull out
	my @start_end = split(';',shift(@ss_location));
	
	foreach my $start_end (@start_end) {
		my @positions = split(',',$start_end);
		$ss_read .= pull_seq($chrom,$positions[0],$positions[1],$index);
	}
	
	### for each particular "splice" of the ss_location (the key of ss_location) that has this sequence, print it out to the map file
	foreach my $ss_splice (keys %{$splice_sites{$ss_location}}) {
		my @keys = sort keys %{$splice_sites{$ss_location}{$ss_splice}};
		### check if sequence is empty
		if ($ss_read eq '') {
			print STDERR "[STDERR]: sequence is empty: $ss_location |", join(',',@keys),"\n";
		} elsif ($ss_read =~ /N/) {
			print STDERR "[STDERR]: sequence has 'N': $ss_location |", join(',',@keys),"\n";
		} else {
			my $display_id = "ss" . $ss_index++;
			print MAP $display_id, "\t", $ss_splice, "\t", join(',',@keys), "\n";
			my $seq_obj = Bio::Seq->new(-seq => $ss_read, -display_id => $display_id, -desc => "", -alphabet => "dna");
			$fasta->write_seq($seq_obj);
		}
	}
}


sub pull_seq {
	my $chrom = shift(@_);
	my $start = shift(@_);
	my $end = shift(@_);
	my $index = shift(@_);

	my $region = $chrom . ":" . $start . "-" . $end;
	my $seq = `samtools faidx $index '$region'`;
	my @split = split('\n',$seq);
	shift(@split);
	return join('',@split);
}
