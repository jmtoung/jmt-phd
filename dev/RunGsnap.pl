#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
umask 007;

my $fastq;
my $fastq2;
my $index;
my $tag;
my $num_threads;
my $part;
my $bam;
my $gsnap_path;
my $samtools_path;
my $B = 4;
my $bufferSize = 1000;
my $snpIndex;
my $spliceIndex;

my $options = GetOptions(
	"fastq=s" => \$fastq, 
	"fastq2=s" => \$fastq2,
	"index=s" => \$index,
	"tag=s" => \$tag,
	"num_threads=i" => \$num_threads,
	"part=s" => \$part,
	"bam=s" => \$bam,
	"gsnap_path=s" => \$gsnap_path,
	"samtools_path=s" => \$samtools_path,
	"B=i" => \$B,
	"bufferSize=i" => \$bufferSize,
	"snpIndex=s" => \$snpIndex,
	"spliceIndex=s" =>\ $spliceIndex
);

$|++;

defined $fastq or die "[STDERR]: fastq not defined\n";
defined $index or die "[STDERR]: index not defined\n";
defined $tag or die "[STDERR]: tag not defined\n" if !defined $bam;
defined $num_threads or die "[STDERR]: threads not defined\n"; 

### split fastq file name ######################################################
my ($fastq_name, $fastq_dir, $fastq_ext) = fileparse($fastq,'\.fastq.gz');

my @part = split('\/',$part) if $part;

unless (defined $bam) {
	$bam = $fastq_dir . $fastq_name;
	$bam =~ s/\/fastq\//\/aln\//;
	$bam =~ s/\/oldhome//; ### added on 4/23 for ir-er time course and onward, get rid of oldhome
	
###	$bam = $bam . ".restart";
	
	unless (-d $bam) { 
		!system("mkdir -p $bam") or die "[STDERR]: $bam doesn't exist\n";
	}
	$bam .= "/gsnap.$fastq_name" . ".$tag";
	if ($part) {
		my $part_tag = 0 x (4 - length($part[0])) . $part[0];
		$bam .= ".$part_tag";
	}
	$bam .= ".bam";
}

### run gsnap ##################################################################
### -B 4 stands for allocation of offsets and positions and genome into memory
print "---running the following command---\n";
print "---start at ", get_timestamp(), "---\n";
my $command = "gsnap -d $index --gunzip -B $B ";
$command .= "-v $snpIndex " if defined $snpIndex;
$command .= "-s $spliceIndex " if defined $spliceIndex;
$command .= "-t $num_threads -N 1 -M 1 -n 10 -Q -O -A sam --input-buffer-size=$bufferSize ";
$command = $gsnap_path . "/" . $command if defined $gsnap_path;
$command .= "--part=$part " if defined $part;
###$command .= " $fastq | samtools view -bS - > $bam"; ### changed on 3/30/2012 to allow for paired-end mapping
$command .= "$fastq ";
$command .= "$fastq2 " if defined $fastq2;
$command .= "| "; ### added samtools path on 3/30
$command .= $samtools_path . "/" if $samtools_path;
$command .= "samtools view -bS - > $bam";

print $command, "\n";
!system($command) or die "[STDERR]: can't run $command: $!\n";
###!system($command) or die "[STDERR]: error!\n";
#system("gsnap -d $index --gunzip -B 4 -v CEU-snp-index -t $num_threads -s refGene_wgEncodeGencodeManualV3 -N 1 -M 1 -n 10 -Q -O -A sam $fastq | samtools view -bS - > $bam");
print "---end at ", get_timestamp(), "---\n";
print "-----------------------------------\n";

print "completed_running:\t$bam\n";

### get time ###################################################################
sub get_timestamp {
        my $TIME = localtime;
        my @TIME = split('\s+',$TIME);
        
        ### get rid of day of week
        shift(@TIME);

        ### convert from Jan to January
        $TIME[0] = full_month($TIME[0]);
        
        return join(" ",@TIME);
}

sub full_month {
        my $MONTH = shift;
        return 'January' if $MONTH eq 'Jan';
        return 'February' if $MONTH eq 'Feb';
        return 'March' if $MONTH eq 'Mar';
        return 'April' if $MONTH eq 'Apr';
        return 'May' if $MONTH eq 'May';
        return 'June' if $MONTH eq 'Jun';
        return 'July' if $MONTH eq 'Jul';
        return 'August' if $MONTH eq 'Aug';
        return 'September' if $MONTH eq 'Sep';
        return 'October' if $MONTH eq 'Oct';
        return 'November' if $MONTH eq 'Nov';
        return 'December' if $MONTH eq 'Dec';
        return $MONTH;
}

