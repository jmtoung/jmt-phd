#!/usr/bin/perl

use lib '/ifs/h/toung/dev', '/home/jmtoung/Lab/dev';

use Database;
use File::Basename;
use Getopt::Long;

my $HOME = "/ifs/h/toung";

my $sites = "/ifs/h/toung/work/2011-09-11_MergeNonRefSites/bowtie_groseq/MergeNonRefSites_bowtie_groseq_rm_poly_min_rdd_reads_5.txt";
my $num_per_file = 10000;
my @bam = qw(
/ifs/h/toung/aln/gro-seq/GM12004PROseqSmLg_trimadapter_trimlowqual/bowtie.GM12004PROseqSmLg_trimadapter_trimlowqual.n2e120.unique.bam
/ifs/h/toung/aln/gro-seq/GM12004PROseqLg_trimadapter_trimlowqual/bowtie.GM12004PROseqLg_trimadapter_trimlowqual.n2e120.unique.bam
/ifs/h/toung/aln/gro-seq/GM12004GRO_FlowCell634MK_N6_trimadapter_trimlowqual/bowtie.GM12004GRO_FlowCell634MK_N6_trimadapter_trimlowqual.n2e120.unique.bam
/ifs/h/toung/aln/adar/GM12004_baseline_trimlowqual/bowtie.GM12004_baseline_trimlowqual.n2e120.unique.bam
);
my $tag = "bowtie_groseq_rm_poly_min_rdd_reads_5";
my $alnDB = "/ifs/h/toung/aln/alnDB";
my $fastqDB = "/ifs/h/toung/fastq/fastqDB";
my $samplesDB = "/ifs/h/toung/database/samplesDB";
my %index;
$index{'F'} = "/ifs/h/toung/database/human_b36_female_ebv.fa";
$index{'M'} = "/ifs/h/toung/database/human_b36_male_ebv.fa";
my $strand_specific = 1;

### divide up sites file into regions #########################################
my $start = 1;
my $end;
my $sites_count = get_file_lines($sites);
my %regions;
my $region_index = 1;
while($start < $sites_count) {
	$end = $start + $num_per_file - 1;
	$regions{$region_index}{'start'} = $start;
	$regions{$region_index}{'end'} = $end;
	$start = $end + 1;
	$region_index++;
}
################################################################################

my @REGIONS;
my @BAM;
my @INDEX;

### foreach bam file... ########################################################
foreach my $bam (@bam) {
	### check that $BAM ends in '.bam' and is an absolute file name ########
	my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
	$bam_ext eq '.bam' or die "[STDERR]: '$BAM' does not end in '.bam'\n";
	substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$BAM' is not absolute file name\n";

	-e $bam or die "[STDERR]: $bam doesn't exist\n";

	### determine sex of individual from bam ###############################
	my $bam_rm_unique = $bam;
	$bam_rm_unique =~ s/\.unique//;
	my $fastqID = Database->new($alnDB)->lookup(1,$bam_rm_unique,2);
	defined $fastqID or die "[STDERR]: fastqID not defined for $bam_rm_unique\n";
	my $sampleID = Database->new($fastqDB)->lookup(0,$fastqID,3);
	defined $sampleID or die "[STDERR]: sampleID not defined for $fastqID\n";
	my $sex = Database->new($samplesDB)->lookup(0,$sampleID,1);
	defined $sex or die "[STDERR]: sex not defined for $sampleID\n";

	foreach my $region_index (sort {$a<=>$b} keys %regions) {
		my $region = $region_index . ":" . $regions{$region_index}{'start'} . "-" . $regions{$region_index}{'end'};
		push(@REGIONS,$region);
		push(@BAM,$bam);
		push(@INDEX,$index{$sex});
	}
	
}

my $region = join("\n",@REGIONS);
my $bam = join("\n",@BAM);
my $index = join("\n",@INDEX);

my $num_files = scalar(@REGIONS);

print STDOUT <<"END";
#\$ -cwd ### use current directory
#\$ -S /usr/bin/perl ### program to execute script
#\$ -M toung\@mail.med.upenn.edu ### email address
##\$ -m ea ### mail is to be sent at abort and end time
#\$ -j y ### combine stdout and stderr
#\$ -t 1-$num_files ### array job #
#\$ -V ### use current environment variables
##\$ -pe DJ 12 ### parallel threads
#\$ -l mem_free=2G ### request memory

use strict;
my \$ID = \$ENV{SGE_TASK_ID} - 1;

my \$sites = "$sites";

my \@region = qw(
$region
);

my \@bam = qw(
$bam
);

my \@index = qw(
$index
);

print "perl $HOME/dev/ExtractRDDInfo.pl --sites \$sites --region \$region[\$ID] --bam \$bam[\$ID] --index \$index[\$ID] --tag $tag --alnDB $alnDB --strand_specific $strand_specific\\n";

system("perl $HOME/dev/ExtractRDDInfo.pl --sites \$sites --region \$region[\$ID] --bam \$bam[\$ID] --index \$index[\$ID] --tag $tag --alnDB $alnDB --strand_specific $strand_specific");

print "completed\\t", "\$sites", "\\t", "\$region[\$ID]", "\\t", "\$bam[\$ID]", "\\t", "\$index[\$ID]", "\\t", \$ID, "\\n";

END



sub get_file_lines {
	my $file = shift;

	my $wc = `wc -l $file`;
	chomp $wc;
	my @wc = split('\s+',$wc);
	return $wc[0];
}
