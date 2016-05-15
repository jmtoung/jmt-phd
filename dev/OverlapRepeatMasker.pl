#!/usr/bin/perl -w

use strict;
use lib '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use File::Basename;
use Getopt::Long;
use RegularList;
use ComplementBase;

my $rdd_file;
my $rmsk_file;

my $result = GetOptions(
	"rdd_file=s" => \$rdd_file,
	"rmsk_file=s" => \$rmsk_file
);

my ($rdd_name, $rdd_dir, $rdd_ext) = fileparse($rdd_file,'\.rdd$') or die "[STDERR]: can't fileparse $rdd_file: $!\n";
$rdd_ext eq '.rdd' or die "[STDERR]: rdd file must end in .rdd\n";
open(RDD_FILE,$rdd_file) or die "can't open $rdd_file: $!\n";
my @rdd_file = split('\/',$rdd_file);
my @rdd_name = split('\.',$rdd_file[-1]);
my $chrom = $rdd_name[-3];

$rmsk_file =~ s/XXX/$chrom/;
-e $rmsk_file or die "[STDERR]: $rmsk_file doesn't exist\n";
my $RMSK = RegularList->new($rmsk_file);

my $output_file = $rdd_file . "_rmsk-temp";
open(OUTPUT,">".$output_file) or die "can't open $output_file: $!\n";
select(OUTPUT);

while(<RDD_FILE>) {
	chomp;
	my @split = split('\t');

	while($RMSK->has_next && $RMSK->get_column(2) < $split[1] && $RMSK->get_column(3) < $split[1]) { $RMSK->update_next_line; }

	### If RDD is within the boundaries of the current gene record	
	if($RMSK->has_next && $split[1] >= $RMSK->get_column(2) && $split[1] <= $RMSK->get_column(3)) {
		
		$split[19] = 1;
	} else {

		$split[19] = 0;
	}
	print join("\t",@split), "\n";
}

system("mv $output_file $rdd_file");
