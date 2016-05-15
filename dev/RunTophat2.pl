#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
umask 007;

my $fastq;
my $fastq2;
my $index;
###my $bam;
my $outputDir;
my $tophat_path;
###my $samtools_path;
my $num_mismatches;
my $num_edit_dist;
my $library_type = "fr-unstranded";
my $transcriptome_index;
my $mateInnerDist;
my $report_secondary = 1;
my $num_threads = 1;

my $options = GetOptions(
	"fastq=s" => \$fastq, 
	"fastq2=s" => \$fastq2,
	"index=s" => \$index,
	"outputDir=s" => \$outputDir,
	"tophat_path=s" => \$tophat_path,
	"num_mismatches=i" => \$num_mismatches,
	"num_edit_dist=i" => \$num_edit_dist,
	"library_type=s" => \$library_type,
	"transcriptome_index=s" => \$transcriptome_index,
	"mateInnerDist=i" => \$mateInnerDist,
	"report_secondary=s" => \$report_secondary,
	"num_threads=i" => \$num_threads
);

$|++;

(defined $fastq && -e $fastq) or die "[STDERR]: fastq not defined\n";
(defined $index) or die "[STDERR]: index not defined\n";
(-d $outputDir && defined $outputDir) or die "[STDERR]: outputDir not defined\n";
(defined $mateInnerDist && $mateInnerDist =~ /^[0-9]+$/) or die "[STDERR]: mateInnerDist not defined\n";
(defined $num_threads && $num_threads =~ /^[0-9]+$/) or die "[STDERR]: threads not defined\n"; 

### run tophat ##################################################################
print "---running the following command---\n";
print "---start at ", get_timestamp(), "---\n";
my $command = "tophat2 -o $outputDir ";
$command .= "-N $num_mismatches " if defined $num_mismatches;
$command .= "--read-edit-dist $num_edit_dist " if defined $num_edit_dist;
$command .= "-p $num_threads " if defined $num_threads;
$command = $tophat_path . "/" . $command if defined $tophat_path;
$command .= "--library-type $library_type " if defined $library_type;
$command .= "--transcriptome-index=$transcriptome_index " if defined $transcriptome_index;
$command .= "-r $mateInnerDist " if defined $mateInnerDist;
$command .= "--report-secondary-alignments " if defined $report_secondary;
$command .= "$index ";
$command .= "$fastq ";
$command .= "$fastq2 " if $fastq2;
###$command .= " | "; ### added samtools path on 3/30
###$command .= $samtools_path . "/" if $samtools_path;
###$command .= "samtools view -bS - > $bam";

print $command, "\n";
!system($command) or die "[STDERR]: can't run $command: $!\n";

print "---end at ", get_timestamp(), "---\n";
print "-----------------------------------\n";

print "completed_running:\t$fastq\n";

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

