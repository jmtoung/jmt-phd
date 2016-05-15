#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
umask 007;

my $fastq;
my $fastq2;
my $index;
my $outputDir;
my $num_chunks = 1;
my $tag;
my $rum_path;
my $qsub = 1;

my $options = GetOptions(
	"fastq=s" => \$fastq, 
	"fastq2=s" => \$fastq2,
	"index=s" => \$index,
	"outputDir=s" => \$outputDir,
	"num_chunks=i" => \$num_chunks,
	"tag=s" => \$tag,
	"rum_path=s" => \$rum_path,
	"qsub=i" => \$qsub
);

$|++;

(defined $fastq && -e $fastq) or die "[STDERR]: fastq not defined\n";
(defined $index) or die "[STDERR]: index not defined\n";
(defined $outputDir && -d $outputDir) or die "[STDERR]: outputDir not defined\n";
(defined $num_chunks && $num_chunks =~ /^[0-9]+$/) or die "[STDERR]: threads not defined\n"; 
(defined $tag) or die "[STDERR]: tag not defined\n";

###my $config = $rum_path . "/conf/rum.config_" . $organism;
###-e $config or die "[STDERR]: config file $config doesn't exist: $!\n";

### run RUM ##################################################################
print "---running the following command---\n";
print "---start at ", get_timestamp(), "---\n";
###my @fastq = ($fastq);
###push(@fastq,$fastq2) if defined $fastq2;
###my $files = join(",,,",@fastq);

my $command = "$rum_path/rum_runner align --output $outputDir --name $tag --index $index --chunks $num_chunks --preserve-names";
####my $command = "perl $rum_path/bin/RUM_runner.pl $config $files $outputDir $num_chunks $tag -preserve_names";
$command .= " --qsub" if $qsub;
$command .= " $fastq";
$command .= " $fastq2" if defined $fastq2;

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

