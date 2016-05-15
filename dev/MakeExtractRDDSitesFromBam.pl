#!/usr/bin/perl

use lib '/ifs/h/toung/dev', '/home/jmtoung/Lab/dev';

use Database;
use File::Basename;
use Getopt::Long;

my $HOME = "/ifs/h/toung";

my $bam;
my $index;
my $tag;
my $exclude;
my $min_rdd_level;
my $min_total_reads;
my $unique_only;
my $adapter_only;
my $strand_specific;
my $alnDB;
my $all_chrom;
my $sge_directory;

my $options = GetOptions ("bam=s" => \$bam, "index=s" => \$index, "tag=s" => \$tag,
			"exclude=s" => \$exclude, "min_rdd_level=i" => \$min_rdd_level,
			"min_total_reads=i" => \$min_total_reads, "unique_only=i" => \$unique_only,
			"adapter_only=i" => \$adapter_only, "strand_specific=i" => \$strand_specific,
			"alnDB=s" => \$alnDB, "all_chrom=i" =>\$all_chrom, "sge_directory=s" =>\$sge_directory);

### check that $BAM ends in '.bam' and is an absolute file name ################
my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
$bam_ext eq '.bam' or die "[STDERR]: '$BAM' does not end in '.bam'\n";
substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$BAM' is not absolute file name\n";

### get chromosomes in $bam ####################################################
my %CHROM;
my $bam_bai = "samtools idxstats $bam |";
open(BAM_BAI,$bam_bai) or die "[STDERR]: can't fork '$bam_bai': $!\n";
while(<BAM_BAI>) {
	chomp;
	my @split = split('\t');
	$CHROM{$split[0]}++;
}
################################################################################

### get chromosomes in $index ##################################################
my $index_fai = $index . ".fai";
open(INDEX_FAI,$index_fai) or die "[STDERR]: can't open fai_index: $index_fai: $!\n";
while(<INDEX_FAI>) {
	chomp;
	my @split = split('\t');
	$CHROM{$split[0]}++;
}
################################################################################

print ">>> Looping through chromosomes <<<\n";
foreach my $CHROM (keys %CHROM) {

	next unless $CHROM{$CHROM} == 2;
	unless($all_chrom) { next unless $CHROM =~ /^(chr)?([0-9XY]+$|^MT$)/; }
	print "\t...working on chromosome '$CHROM'...\n";

	### get EXCLUDE list ready #############################################
	my @exclude_new;
	foreach my $EXCLUDE (split(',',$exclude)) {
		my $exclude_new = $EXCLUDE;
		$exclude_new =~ s/XXX/$CHROM/;
		if(-e $exclude_new) { push(@exclude_new,$exclude_new); }
		else { print "\t\t[WARNING]: $exclude_new doesn't exist\n"; }
	}
	my $exclude_new = join(',',@exclude_new);

	### print out files ####################################################
	my $RUN_FILE = $sge_directory . "/sge_ExtractRDDSitesFromBam_${bam_name}_${tag}_$CHROM.sh"; 

	open(RUN_FILE,">$RUN_FILE") or die "[STDERR]: can't open $RUN_FILE: $!\n";

print RUN_FILE <<"END";
#\$ -cwd ### use current directory
#\$ -S /usr/bin/perl ### program to execute script
#\$ -M toung\@mail.med.upenn.edu ### email address
##\$ -m ea ### mail is to be sent at abort and end time
#\$ -j y ### combine stdout and stderr
#\$ -t 1-1 ### array job #
#\$ -V ### use current environment variables
##\$ -pe DJ 12 ### parallel threads
#\$ -l mem_free=4G ### request memory

my \$ID = \$ENV{SGE_TASK_ID} - 1;

perl $HOME/dev/ExtractRDDSitesFromBam.pl --bam $bam --index $index --chrom $chrom --tag $tag --exclude $exclude --min_rdd_level $min_rdd_level --min_total_read $min_total_reads --unique_only $unique_only --adapter_only $adapter_only --strand_specific $strand_specific

END
	close(RUN_RILE);
	}
}
