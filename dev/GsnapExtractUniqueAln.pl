#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Bio::DB::Sam;
use File::Basename;
use Statistics::Descriptive;
use RegularList;
use PileupData;
use ComplementBase;

my $HOME = "/ifs/h/toung";
my $bam;
my $min_mapqual;

my $options = GetOptions(
	"home=s" => \$HOME,
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
### make output
my $output_file = $bam_dir . "$bam_name.unique$bam_ext";
################################################################################

system("perl $HOME/dev/GsnapExtractUniqueAlnWrapper.pl --bam $bam --min_mapqual $min_mapqual | samtools view -bS - > $output_file");

################################################################################################################
### this was copied over from the ~/jobs/RunGsnap_def/step5_SamtoolsSortIndex/4/sge_samtools-index.pl script ###
################################################################################################################

($bam_name,$bam_dir,$bam_ext) = fileparse($output_file,'\.bam');

my $bam_bai = $output_file . ".bai";
my $bam_sort = $bam_dir . $bam_name . "_sort" . $bam_ext;
my $sort_name = $bam_name . "_sort";

print "bamfile=>\t", $output_file, "\n";
print "bam_bai_file=>\t", $bam_bai, "\n";
print "bam_sort_file=>\t", $bam_sort, "\n";
print "sort_name=>\t", $sort_name, "\n";

chdir($bam_dir);

system("samtools sort $output_file $sort_name");
system("mv $bam_sort $output_file");
system("samtools index $output_file");

print "sorting_success_for\t$output_file\n";
