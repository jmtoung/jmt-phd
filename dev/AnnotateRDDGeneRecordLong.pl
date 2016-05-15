#!/usr/bin/perl -w

use strict;
use lib '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use File::Basename;
use Getopt::Long;
use RegularList;
use ComplementBase;

my $rdd_file = "/home/jmtoung/Lab/aln/adar/GM12750_ADAR_KD_trimlowqual/bowtie.GM12750_ADAR_KD_trimlowqual.n2e120.unique.minqual20-minlvl10-minct10.rdd/bowtie.GM12750_ADAR_KD_trimlowqual.n2e120.unique.minqual20-minlvl10-minct10.Y:1-57772954.rdd";
my $annotation_file = "/home/jmtoung/Lab/database/annotation_files/gencode.v3_refseq_human_b36/gencode.v3_refseq_human_b36.XXX.txt";
my $num_col = 19;

my $result = GetOptions("rdd_file" => \$rdd_file, "annotation_file=s" => \$annotation_file, "num_col=i" => \$num_col);

my ($rdd_name, $rdd_dir, $rdd_ext) = fileparse($rdd_file,'\.rdd') or die "[STDERR]: can't fileparse $rdd_file: $!\n";
$rdd_ext eq '.rdd' or die "[STDERR]: rdd file must end in .rdd\n";
open(RDD_FILE,$rdd_file) or die "can't open $rdd_file: $!\n";
my @rdd_file = split('\/',$rdd_file);
my @rdd_name = split('\.',$rdd_file[-1]);
my $region = $rdd_name[-2];
my @chrom = split(':',$region);

$annotation_file =~ s/XXX/$chrom[0]/;
-e $annotation_file or die "[STDERR]: $annotation_file doesn't exist\n";
my $ANNOTATION = RegularList->new($annotation_file);

my $output_file = $rdd_dir . $rdd_name . ".annotate" . $rdd_ext;
open(OUTPUT,">".$output_file) or die "can't open $output_file: $!\n";
select(OUTPUT);

while(<RDD_FILE>) {
	chomp;
	my @split = split('\t');
	
	while(scalar(@split) != $num_col) {
		scalar(@split) > $num_col or die "[STDERR]: fewer col than expected\n";
		pop(@split);
	}

	while($ANNOTATION->has_next && $ANNOTATION->get_column(1) < $split[1] && $ANNOTATION->get_column(2) < $split[1]) { $ANNOTATION->update_next_line; }

	### If RDD is within the boundaries of the current gene record	
	if($split[1] >= $ANNOTATION->get_column(1) && $split[1] <= $ANNOTATION->get_column(2)) {
		my @GENE_RECORDS = split(';',$ANNOTATION->get_column(3));
		my %ANNOTATION;
		foreach my $GENE_RECORD (@GENE_RECORDS) {
			my ($ENSG_ID, $STRAND, $TYPE) = split('\|',$GENE_RECORD);
			foreach my $type (split(',',$TYPE)) {
				$ANNOTATION{$STRAND}{$ENSG_ID}{$type}++;
			}
		}

		my @STRANDS = keys %ANNOTATION;

		### If all gene records are on one and only one strand
		if (scalar(@STRANDS) == 1) {
			if ($STRANDS[0] eq '+') {
				if ($split[2] ne '-') {
					$split[2] = '+';
					foreach my $ENSG_ID (keys %{$ANNOTATION{$STRANDS[0]}}) {
						foreach my $TYPE (keys %{$ANNOTATION{$STRANDS[0]}{$ENSG_ID}}) {
							print join("\t",@split), "\t", $ENSG_ID, "\t", $TYPE, "\n";
						}
					}
				} else {
					print STDERR "[STDERR]: discordant strand: $_\n";
					$split[2] = '?';
					print join("\t",@split), "\t\t\n";
				}
			} else {
				if ($split[2] ne '+') {
					$split[2] = '-';
					$split[3] = complement_base($split[3]);
					$split[4] = complement_base($split[4]);
					foreach my $ENSG_ID (keys %{$ANNOTATION{$STRANDS[0]}}) {
						foreach my $TYPE (keys %{$ANNOTATION{$STRANDS[0]}{$ENSG_ID}}) {
							print join("\t",@split), "\t", $ENSG_ID, "\t", $TYPE, "\n";
						}
					}
				} else {
					print STDERR "[STDERR]: discordant strand: $_\n";
					$split[2] = '?';
					print join("\t",@split), "\t\t\n";
				}
			}
		### If all gene records are on two strands
		} elsif (scalar(@STRANDS) == 2) {
			foreach my $STRAND (@STRANDS) {
				if (($split[2] eq '+' && $STRAND ne '+') || ($split[2] eq '-' && $STRAND ne '-')) { 
					print STDERR "[STDERR]: discordant strand: $_\n"; 
					$split[2] = '?';
					print join("\t",@split), "\t\t\n";
					next;					
				}
				$split[2] = $STRAND;
				foreach my $ENSG_ID (keys %{$ANNOTATION{$STRAND}}) {
					foreach my $TYPE (keys %{$ANNOTATION{$STRAND}{$ENSG_ID}}) {
						print join("\t",@split), "\t", $ENSG_ID, "\t", $TYPE, "\n";
					}			
				}
			}
		}
	### If RDD is not within current gene record (meaning it's not in list at all -- intergenic?)
	} else {
		print join("\t",@split), "\t\tintergenic\n", 
	}
}
