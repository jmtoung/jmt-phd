#$ -cwd ### use current directory
#$ -S /usr/bin/perl ### program to execute script
#$ -M toung@mail.med.upenn.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-3 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
##$ -l mem_free=6G ### request memory

my $ID = $ENV{SGE_TASK_ID} - 1;

use lib '/ifs/h/toung/dev', '/home/jmtoung/Lab/';
use Database;
use strict;

### CHANGE THESE VARIABLES
my $HOME = '/ifs/h/toung';
my $SAMPLE = '12750';
my $SAMPLESDB_FILE = $HOME . "/database/samplesDB.txt";
my $INDEX_FEMALE = $HOME . "/database/human_b36_female_ebv.fa";
my $INDEX_MALE = $HOME . "/database/human_b36_male_ebv.fa";
my $EXPERIMENT = "adar";
my @FASTQ_FILE = qw(
GM12750_ADAR_KD_trimlowqual
GM12750_NT_1_1_trimlowqual
GM12750_NT_1_2_trimlowqual
);
my $FASTQ_FILE = $FASTQ_FILE[$ID];
my $BAM_FILE = $HOME . "/aln/$EXPERIMENT/$FASTQ_FILE/bowtie.$FASTQ_FILE.n2e120.bam";
my $EXCLUDE = "/ifs/h/toung/rdd/1000_genomes/pilot_data/release/2010_07/1000G-INDEL.XXX.txt,/ifs/h/toung/rdd/1000_genomes/pilot_data/release/2010_07/1000G-SNP.XXX.txt,/ifs/h/toung/rdd/affymetrix/affymetrix_all.XXX.txt,/ifs/h/toung/rdd/hapmap/genotypes/2010-08_phaseII+III/forward/hapmap-release28.XXX.txt";
my $TAG = "unique-10pct";
my $MIN_RDD_LEVEL = .1;
my $MIN_TOTAL_READS = 0;
my $UNIQUE_ONLY = 1;
my $ADAPTER_ONLY = 0;
my $STRAND_SPECIFIC = 0;
my $ALNDB = $HOME . "/aln/alnDB/alnDB.txt";
my $ALL_CHROM = 0;
my $SGE_DIRECTORY = $HOME . "/aln/$EXPERIMENT";
$|++;

### LOOK UP SEX OF SAMPLE ######################################################
my $SAMPLESDB = Database->new($SAMPLESDB_FILE);
my $SEX = $SAMPLESDB->lookup_val(0,$SAMPLE,1);
defined $SEX or die "[STDERR]: $SAMPLE doesn't exist in $SAMPLESDB_FILE\n";
################################################################################

### GET INDEX ##################################################################
my $INDEX;
if ($SEX eq 'F') { $INDEX = $INDEX_FEMALE }
else { $INDEX = $INDEX_MALE; }
################################################################################

print "perl $HOME/dev/MakeExtractRDDSitesFromBam.pl --bam $BAM_FILE --index $INDEX --tag $TAG --exclude $EXCLUDE";

system("perl $HOME/dev/MakeExtractRDDSitesFromBam.pl --bam $BAM_FILE --index $INDEX --tag $TAG --exclude $EXCLUDE --min_rdd_level $MIN_RDD_LEVEL --min_total_reads $MIN_TOTAL_READS --unique_only $UNIQUE_ONLY --adapter_only $ADAPTER_ONLY --strand_specific $STRAND_SPECIFIC --alnDB $ALNDB --all_chrom $ALL_CHROM --sge_directory $SGE_DIRECTORY");
