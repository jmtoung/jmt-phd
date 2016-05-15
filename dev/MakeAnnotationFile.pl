#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];

### gtf files, separated by commas
my $gtf = "/home/jmtoung/Lab/rdd/1000_genomes/technical/retired_reference/gencode.v3c.annotation.NCBI36_chr_removed.gtf,/home/jmtoung/Lab/database/refseqgenes_refgene_03072011_chr_removed.gtf";
my $output_folder = "/home/jmtoung/Lab/database/annotation_files/gencode.v3_refseq_human_b36";
my $tag = "gencode.v3_refseq_human_b36";

my $result = GetOptions("gtf=s" => \$gtf,"output_folder=s" => \$output_folder, "tag=s" => \$tag);

#######################################################################################
### split gtf
my @gtf = split(',',$gtf);
### make directory
unless(-d $output_folder) { mkdir $output_folder; }
chdir($output_folder);
#######################################################################################

###  GROUP RECORDS IN GTF FILE BY THE "transcript_id" in %gtf hash
my %gtf;
foreach my $GTF (@gtf) {
	
	open(GTF,$GTF) or die "[STDERR:] error opening gtf file $GTF: $!\n";
	### this is the variable we will use to group features of the same "transcript". should be "transcript_id".
	my $index_var = "transcript_id";
	while(<GTF>) {
		next if (/^#/);
		chomp;
		### split the line
		my @split = split('\t');

		### we are going to split column 9 (last col), the attributes column and store the keys and values in a hash %split2
		my %split2;
		foreach (split(';',pop(@split))) {
			if (/\s*(.*)\s"(.*)"/) {
				$split2{$1} = $2;
			} elsif (/level\s([0-9])+/) {
				$split2{"level"} = $1;
			}
		}
		
		### store variables into hash
		my $index = $split2{$index_var};
		### chrom->strand->transcript_id->type->start->end
		### make sure that if type is 'gene' or 'transcript', only one exists
		if ($split[2] eq 'gene' || $split[2] eq 'transcript') {
			if (exists $gtf{$split[0]}{$split[6]}{$index}{$split[2]}) { print STDERR "error: gene/transcript has two entries: ", $_, "\n"; }
			else { $gtf{$split[0]}{$split[6]}{$index}{$split[2]}{$split[3]} = $split[4]; }	
		### for anything else, make sure that there aren't two ends for a start
		} else {
			if (exists $gtf{$split[0]}{$split[6]}{$index}{$split[2]}{$split[3]}) { print STDERR "error: same type two starts: ", $_, "\n"; }
			else { $gtf{$split[0]}{$split[6]}{$index}{$split[2]}{$split[3]} = $split[4]; }
		}
	}
	close(GTF);
}

### loop through each chromosome and tally up/compress all records 
foreach my $chrom (keys %gtf) {
	
	my $file = $tag . '.' . $chrom . '.txt';
	open(OUTPUT,'>'.$file) or die "[STDERR]: can't open $file: $!\n";
	select(OUTPUT);

	### %GTF is a hash that has the start of the "record" as a key and an array with the end, type, and strand
	### store the starts of ever record on this chromosome here.
	my %GTF = ();

	### LOOP THROUGH ALL RECORDS ###########################################
	foreach my $strand (keys %{$gtf{$chrom}}) {
		foreach my $transcript (keys %{$gtf{$chrom}{$strand}}) {
			
			### this is a hash that keeps track of exons so we can determine introns
			my %structure;
			my (@STARTS,@ENDS); ### use this to store the start and end for the entire gene/transcript. 
			
			### GENE / TRANSCRIPT
			foreach my $type ('gene','transcript') {
				if (exists $gtf{$chrom}{$strand}{$transcript}{$type}) {
					### use the each construct if we don't care about getting entries back in order. in this case we don't since there should only be one entry
					while(my($start,$end) = each (%{$gtf{$chrom}{$strand}{$transcript}{$type}})) {
						push(@{$GTF{$start}},[$end,$transcript,$type,$strand]);
						push(@STARTS,$start);
						push(@ENDS,$end);
					}
				}
			}
			
			### EXON
			### (add splice sites here too)
			if (exists $gtf{$chrom}{$strand}{$transcript}{'exon'}) {
				### get @starts so that we can get exons in order
				my @starts = sort {$a <=> $b} keys %{$gtf{$chrom}{$strand}{$transcript}{'exon'}};
				for (my $i = 0; $i <= $#starts; $i++) {
					### get the end of the start
					my $end = $gtf{$chrom}{$strand}{$transcript}{'exon'}{$starts[$i]};
					push(@{$GTF{$starts[$i]}},[$end,$transcript,'exon',$strand]);
					push(@STARTS,$starts[$i]);
					push(@ENDS,$end);
					### add exon to the structure hash (for determining introns later)
					structure_tally(\%structure,$starts[$i],$end);
				}
				
				### determine splice sites
				for (my $i = 0; $i < $#starts; $i++) {
					### start and end of 3' end of splice site
					
					my $start = $gtf{$chrom}{$strand}{$transcript}{'exon'}{$starts[$i]} + 1;
					my $end = $start + 1;
					push(@{$GTF{$start}},[$end,$transcript,'ss',$strand]);
					### start and end of 5' end of splice site
					$start = $starts[$i + 1] - 2;
					$end = $start + 1;
					push(@{$GTF{$start}},[$end,$transcript,'ss',$strand]);
					### don't add this to structure, since it's intronic
				}
			}
			
			### CDS
			my $min_CDS_start; ### save the min start position of the CDS to determine UTR directionality
			if (exists $gtf{$chrom}{$strand}{$transcript}{'CDS'}) {
			my @starts = sort {$a <=> $b} keys %{$gtf{$chrom}{$strand}{$transcript}{'CDS'}};
				$min_CDS_start = $starts[0];
				for (my $i = 0; $i <= $#starts; $i++) {
					my $end = $gtf{$chrom}{$strand}{$transcript}{'CDS'}{$starts[$i]};
					push(@{$GTF{$starts[$i]}},[$end,$transcript,'CDS',$strand]);
					push(@STARTS,$starts[$i]);
					push(@ENDS,$end);					
					structure_tally(\%structure,$starts[$i],$end);
				}
			}
			
			### UTR
			if (exists $gtf{$chrom}{$strand}{$transcript}{'UTR'}) {
				my @starts = sort {$a <=> $b} keys %{$gtf{$chrom}{$strand}{$transcript}{'UTR'}};
				for (my $i = 0; $i <= $#starts; $i++) {
					my $end = $gtf{$chrom}{$strand}{$transcript}{'UTR'}{$starts[$i]};
					push(@STARTS,$starts[$i]);
					push(@ENDS,$end);
					if ($min_CDS_start) { ### use the minimum CDS start as an anchor to determine 3' versus 5' UTR
						if (($strand eq '+' && $starts[$i] < $min_CDS_start) || ($strand eq '-' && $starts[$i] > $min_CDS_start)) {
							push(@{$GTF{$starts[$i]}},[$end,$transcript,'5UTR',$strand]);
						} else {
							push(@{$GTF{$starts[$i]}},[$end,$transcript,'3UTR',$strand]);
						}
					} else {
						push(@{$GTF{$starts[$i]}},[$end,$transcript,'UTR',$strand]);	
					}
					structure_tally(\%structure,$starts[$i],$end);
				}				
			
			}
			### START CODON / STOP CODON
			foreach my $type ('start_codon','stop_codon') {
				if (exists $gtf{$chrom}{$strand}{$transcript}{$type}) {
					while(my($start,$end) = each (%{$gtf{$chrom}{$strand}{$transcript}{$type}})) {
						push(@STARTS,$start);
						push(@ENDS,$end);						
						push(@{$GTF{$start}},[$end,$transcript,$type,$strand]);
						structure_tally(\%structure,$start,$end);
					}
				}
			}
			
			### DETERMINE INTRONS (only do this if %structure is not empty)
			my $START = min(@STARTS);
			my $END = max(@ENDS);
			if(scalar(%structure)) {
				my @INTRON_BOUNDARIES; ### keeps track of the starts and end boundaries in order (alternates between starts and ends)
				my $IS_INTRON = ! exists $structure{$START}; ### truth means it's an intron
				push(@INTRON_BOUNDARIES,$START) if $IS_INTRON; ### start intron_boundaries if it's an intron
				my $NEW_IS_INTRON; ### use this to compare subsequent position to previous
				### loop through all positions (except first and last) as defined by 'gene' or 'transcript'
				for (my $i = ($START + 1); $i <= $END; $i++) {
					$NEW_IS_INTRON = ! exists $structure{$i};
					next if ($NEW_IS_INTRON == $IS_INTRON); ### next if there's no change
					if ($IS_INTRON) { push(@INTRON_BOUNDARIES,$i - 1); } ### if we're going from intron (IS_INTRON) to exon (NEW_IS_INTRON), then previous $i was last "intron" position
					else { push(@INTRON_BOUNDARIES,$i); }
					$IS_INTRON = $NEW_IS_INTRON;
				}
				$NEW_IS_INTRON = ! exists $structure{$END};
				push(@INTRON_BOUNDARIES,$END) if $NEW_IS_INTRON; ### this is the base case.
	
				### while array exists, add boundaries to %GTF
				while(@INTRON_BOUNDARIES) {
					my $start = shift(@INTRON_BOUNDARIES);
					my $end = shift(@INTRON_BOUNDARIES);
					push(@{$GTF{$start}},[$end,$transcript,'intron',$strand]);
				}
			}
		}
	}
	########################################################################

	### LOOP THGOUH ALL STARTS AND PRINT ###################################
	my @STARTS = sort {$a <=> $b} keys %GTF;
	my $POS = ($STARTS[0] - 1);
	my %STACK;
	my ($FINAL_OUTPUT,$NEW_FINAL_OUTPUT,$START);
	while(@STARTS || scalar(%STACK)) {
		$POS++;
		
		while(@STARTS && $STARTS[0] == $POS) {
			my $RECORDS = $GTF{shift(@STARTS)};
			foreach my $record (@{$RECORDS}) { 
				$STACK{$record->[0]}{$record->[3]}{$record->[1]}{$record->[2]}++;
			}
		}
		my %TO_PRINT;
		
		foreach my $END (keys %STACK) {
			if ($END >= $POS) { 
				foreach my $STRAND (keys %{$STACK{$END}}) {
					foreach my $GENE (keys %{$STACK{$END}{$STRAND}}) {
						foreach my $TYPE (keys %{$STACK{$END}{$STRAND}{$GENE}}) {
							$TO_PRINT{$STRAND}{$GENE}{$TYPE}++;
						}
					}
				}
			} else { delete $STACK{$END}; }
		}
		my %TO_PRINT_COMPACT;
		foreach my $STRAND (sort keys %TO_PRINT) {
			foreach my $GENE (sort keys %{$TO_PRINT{$STRAND}}) {
				my @TYPES = sort keys %{$TO_PRINT{$STRAND}{$GENE}};
				$TO_PRINT_COMPACT{$GENE.'|'.$STRAND.'|'.join(',',@TYPES)}++;
			}
		}
		my @FINAL_OUTPUT = sort keys %TO_PRINT_COMPACT;
		$NEW_FINAL_OUTPUT = join(';',@FINAL_OUTPUT);
			if ($FINAL_OUTPUT) {
			if ($NEW_FINAL_OUTPUT ne $FINAL_OUTPUT) {
				print $chrom, "\t", $START, "\t", ($POS - 1), "\t", $FINAL_OUTPUT, "\n";
				$FINAL_OUTPUT = $NEW_FINAL_OUTPUT;
				$START = $POS;
			}
		} else {
			$FINAL_OUTPUT = $NEW_FINAL_OUTPUT;
			$START = $POS;
		}	
	}
	close(OUTPUT);
}

sub structure_tally {
	my $hash = shift;
	my $start = shift;
	my $end = shift;
	
	for (my $i = $start; $i <= $end; $i++) { $hash->{$i}++; }
}
