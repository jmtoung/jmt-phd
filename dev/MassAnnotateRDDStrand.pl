#!/usr/bin/perl

use lib '/ifs/h/toung/dev', '/home/jmtoung/Lab/dev';

use Database;
use File::Basename;
use Getopt::Long;
use File::Find;

my $HOME = "/home/jmtoung/Lab";

my @directory = qw(
/home/jmtoung/Lab/rdd/twins/GM14452_s_3_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14381_s_7_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14433_s_4_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14468_s_2_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14569_s_2_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14432_s_3_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14467_s_1_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14568_s_1_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14507_s_6_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14448_s_8_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14521_s_6_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14506_s_5_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14453_s_4_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14447_s_7_sequence_trimlowqual
/home/jmtoung/Lab/rdd/twins/GM14382_s_8_sequence_trimlowqual
);

our $annotation_file = "/ifs/h/toung/database/annotation_files/gencode.v3_refseq_human_b36/gencode.v3_refseq_human_b36.XXX.txt";

find(\&process_file,@directory);

sub process_file {
	next unless /\.rdd/;
	$DB::single = 1;
	system("perl $HOME/dev/AnnotateRDDStrand.pl --rdd_file $File::Find::name --annotation_file $annotation_file");

}
