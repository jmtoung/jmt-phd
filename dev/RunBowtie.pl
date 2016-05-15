#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use lib '/ifs/h/toung/dev','/home/jmtoung/Lab/dev';
use CalculateBamStats;
use Database;
use Timestamp;

my $HOME = '/home/jmtoung/Lab';

my $fastqDB;
my $fastq;
my $index;
my $bam;
my $options;
my $ss_map;
my $max_edit_distance;
my $max_strata;
my $alnDB;

$|++;

my $option = GetOptions(
	"fastqDB=s" => \$fastqDB, "fastq=s" => \$fastq, 
	"index=s" => \$index, "bam=s" => \$bam, "options=s" => \$options,
	"ss_map=s" => \$ss_map, "max_edit_distance=s" => \$max_edit_distance,
	"max_strata=s" => \$max_strata, "alnDB=s" => \$alnDB);

### unsplit options ############################################################
my @options = split('\.',$options);
my $options_split = join(" ",@options);
################################################################################

### look up fastq in fastqDB, get fastqID ######################################
my $FASTQ_DB = Database->new($fastqDB);
my $fastqID = $FASTQ_DB->lookup_val(1,$fastq,0); ### foreign key
defined $fastqID or die "[STDERR]: '$fastq' doesn't exist in $fastqDB\n";
################################################################################

### check that $bam ends in '.bam' and is an absolute file name ################
my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
$bam_ext eq '.bam' or die "[STDERR]: '$bam' does not end in '.bam'\n";
substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$bam' is not absolute file name\n";
chdir($bam_dir) or die "cannot chdir '$bam_dir': $!\n";
################################################################################

#### make temp files to store stderr of bowtie and samtools commands ############
#my ($temp1, $temp2);
#do {
#	my $rand = int(rand(100));
#	$temp1 = $bam_dir . $bam_name . "_" . $rand . ".temp1";
#	$temp2 = $bam_dir . $bam_name . "_" . $rand . ".temp2";
#} until (!-e $temp1 && !-e $temp2);
#################################################################################

#### run bowtie and convert to bam in one step ##################################
#print "*** Running Bowtie ***\n";
#print "bowtie $options_split $index $fastq 2> $temp1 | samtools view -bS - > $bam 2> $temp2\n";
#!system("bowtie $options_split $index $fastq 2> $temp1 | samtools view -bS - > $bam 2> $temp2") or die "[STDERR]: cannot run bowtie: $!\n";
#print "*** Finished Running Bowtie ***\n";
#### capture stderr from running bowtie and converting to bam
#open(TEMP1,$temp1) or die "[STDERR]: cannot open '$temp1': $!\n";
#open(TEMP2,$temp2) or die "[STDERR]: cannot open '$temp2': $!\n";
#################################################################################

#### from the stderr of the bowtie run, count up number of reads exhausted and runtime
#my $num_warning_exhausted = 0;
#my $runtime;
#while(<TEMP1>) {
#	chomp; 
#	print "bowtie: $_\n";
#	$num_warning_exhausted++ if (/Warning: Exhausted best-first chunk memory for read/);
#	if (/Time searching: ([0-9:]+)/) { $runtime = $1; }
#}
#### if "runtime" is not defined, then bowtie didn't run properly
#defined $runtime or die "[STDERR]: bowtie didn't run properly: $!\n";
#################################################################################

#### just for redirection of stderr #############################################
#while(<TEMP2>) { print "samtools view: $_"; }
#################################################################################
my $num_warning_exhausted = 0;
my $runtime = '03:25:20';
#### convert coordinates ########################################################
#print "*** Converting Splice Read Coordinates ***\n";
!system("perl $HOME/dev/ConvertSpliceReadCoordinates.pl --bam $bam --ss_map $ss_map") or die "[STDERR]: can't convert coordinates: $!\n";
print "*** Finished Converting Splice Read Coordinates ***\n";

### clean bowtie
print "*** Cleaning Bowtie ***\n";
!system("perl $HOME/dev/CleanBowtieBam.pl --bam $bam --max_edit_distance $max_edit_distance --max_strata $max_strata") or die "[STDERR]: can't clean bowtie: $!\n";
print "*** Finished Cleaning Bowtie ***\n";

### create bam index
print "*** Making Bam Index ***\n";
!system("samtools index $bam") or die "[STDERR]: can't make index for $bam: $!\n";
print "*** Finished Making Bam Index ***\n";

### Calculate number aligned ###################################################
my ($num_aligned,$hist_alignment) = calculate_bam_stats($bam);
################################################################################

my @entry = ($bam,$fastqID,$index,"bowtie",sort_options($options),$num_aligned,
	$hist_alignment,$num_warning_exhausted,$runtime);

my $ALN_DB = Database->new($alnDB);
$ALN_DB->add_entry(\@entry,1); ### pkey is bam file

sub sort_options {
	my $options = shift;

	my @options = split('\s+',$options);
	my %options;
	my $value;
	while (@options) {
		my $token = pop(@options);
		if ($token =~ /-/) { 
			$options{$token} = $value; 
			$value = undef;
		} else { 
			$value = $token; 
		}
	}

	foreach my $key (keys %options) { 
		if (defined $options{$key}) {		
			push(@options,$key . " " . $options{$key});
		} else {
			push(@options,$key);
		}
	}

	return join(" ",@options);
}

END {
#	if (-e $temp1) { system("rm $temp1"); }
#	if (-e $temp2) { system("rm $temp2"); }
}

