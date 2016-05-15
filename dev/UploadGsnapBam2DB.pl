#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use lib '/ifs/h/toung/dev','/home/jmtoung/Lab/dev';
use CalculateBamStats;
use Database;
use Timestamp;

umask 0007;

my $HOME; ##  = '/ifs/h/toung'; changed on 12/18/2011
my $fastqDB;
my $bam;
my $index;
my $options; # = "-B.4.-v.CEU-snp-index.-s.refGene_wgEncodeGencodeManualV3.-N.1.-M.1.-n.10.-Q.-O.-A.sam";
my $alnDB;

$|++;

my $option = GetOptions(
	"home=s" => \$HOME,
	"fastqDB=s" => \$fastqDB,
	"bam=s" => \$bam,
	"index=s" => \$index,
	"options=s" => \$options,
	"alnDB=s" => \$alnDB
);

### unsplit options ############################################################
my @options = split('\.',$options);
my $options_split = join(" ",@options);
################################################################################

my ($bam_name,$bam_dir,$bam_ext) = fileparse($bam,'\.bam');
my $fastq = $bam_dir;
$fastq =~ s/aln/fastq/;
chop $fastq;
$fastq = $fastq . ".fastq";
print "fastq:\t", $fastq, "\n";

### remove $HOME from $fastq ### added 12/18/2011
my $HOME_REGEX = $HOME;
$HOME_REGEX =~ s/\//\\\//g;
my $fastq_rmhome = $fastq;
$fastq_rmhome =~ s/$HOME_REGEX//g;
$fastq_rmhome =~ s/^\///;

### look up $fastq in fastqDB to get fastqID ###################################
my $fastqID = Database->new($fastqDB)->lookup(1,$fastq_rmhome,0);
defined $fastqID or die "[STDERR]: FASTQFAIL '$fastq' doesn't exist in $fastqDB\n";
################################################################################

### calculate number aligned ###################################################
my ($num_aligned,$hist_alignment) = calculate_bam_stats($bam);
################################################################################

my $num_warning_exhausted = 0;
my $runtime = '-';

### remove $HOME from $bam
my $bam_rmhome = $bam;
$bam_rmhome =~ s/$HOME_REGEX//g;
$bam_rmhome =~ s/^\///;

### add entry to database ######################################################
my @entry = (	
	$bam_rmhome,
	$fastqID,
	$index,
	"gsnap",
	sort_options($options),
	$num_aligned,
	$hist_alignment,
	$num_warning_exhausted,
	$runtime,
	"'" . get_timestamp . "'"
);
print "***Adding Entry to Database ***\n";
print "\n", join("\t",@entry), "\n";
Database->new($alnDB)->add2db(\@entry,1);
################################################################################

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

