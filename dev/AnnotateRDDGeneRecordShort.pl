#!/usr/bin/perl -w

use strict;
use lib '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use File::Basename;
use Getopt::Long;
use RegularList;
use ComplementBase;

my $rdd_file;
my $annotation_file;

my $result = GetOptions(
	"rdd_file=s" => \$rdd_file,
	"annotation_file=s" => \$annotation_file
);

my ($rdd_name, $rdd_dir, $rdd_ext) = fileparse($rdd_file,'\.rdd$') or die "[STDERR]: can't fileparse $rdd_file: $!\n";
$rdd_ext eq '.rdd' or die "[STDERR]: rdd file must end in .rdd\n";
open(RDD_FILE,$rdd_file) or die "can't open $rdd_file: $!\n";
my @rdd_file = split('\/',$rdd_file);
my @rdd_name = split('\.',$rdd_file[-1]);
my $chrom = $rdd_name[-3];

$annotation_file =~ s/XXX/$chrom/;
-e $annotation_file or die "[STDERR]: $annotation_file doesn't exist\n";
my $ANNOTATION = RegularList->new($annotation_file);

my $output_file = $rdd_file . "_temp";
open(OUTPUT,">".$output_file) or die "can't open $output_file: $!\n";
select(OUTPUT);

while(<RDD_FILE>) {
	chomp;
	my @split = split('\t');
	
	while($ANNOTATION->has_next && $ANNOTATION->get_column(1) < $split[1] && $ANNOTATION->get_column(2) < $split[1]) { $ANNOTATION->update_next_line; }

	### If RDD is within the boundaries of the current gene record	
	if($ANNOTATION->has_next && $split[1] >= $ANNOTATION->get_column(1) && $split[1] <= $ANNOTATION->get_column(2)) {
		my @GENE_RECORDS = split(';',$ANNOTATION->get_column(3));
		my %TYPES;
		foreach my $GENE_RECORD (@GENE_RECORDS) {
			my ($ENSG_ID, $STRAND, $TYPE) = split('\|',$GENE_RECORD);
			foreach my $type (split(',',$TYPE)) {
				next if $type eq 'transcript' or $type eq 'gene';
				if ($type eq 'exon' or $type eq 'intron') { $TYPES{'exon/intron'}{$type}++; }
				else { $TYPES{'other'}{$type}++; }
			}
		}

		my @exon_intron = keys %{$TYPES{'exon/intron'}};
		unless(@exon_intron) { $split[21] = '-'; }
		else { $split[20] = join('+',sort @exon_intron); }

		my @others = keys %{$TYPES{'other'}};
		unless(@others) { $split[21] = '-'; } 		
		else { $split[21] = join('+',sort @others); }
		
	### If RDD is not within current gene record (meaning it's not in list at all -- intergenic?)
	} else {
		$split[20] = 'intergenic';
		$split[21] = 'intergenic';
	}
	
	print join("\t",@split), "\n";
}

system("mv $output_file $rdd_file");
