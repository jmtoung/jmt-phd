#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

my $bam;
my $byReadName = 0;
my $replaceOriginal = 0;
my $samtools = "samtools";

my $options = GetOptions(
	"bam=s" => \$bam,
	"byReadName=i" => \$byReadName,
	"replaceOriginal=i" => \$replaceOriginal,
	"samtools=s" => \$samtools
);

(print STDERR "bam:\t$bam\n") && -e $bam or die "[STDERR]: $bam not defined/exists\n";
(print STDERR "byReadName:\t$byReadName\n") && defined $byReadName or die "[STDERR]: byReadName not defined\n";
(print STDERR "replaceOriginal:\t$replaceOriginal\n") && defined $replaceOriginal or die "[STDERR]: replaceOriginal not defined\n";
print STDERR "samtools:\t$samtools\n";

my ($bamName,$bamDir,$bamExt) = fileparse($bam,'\.bam');
$bamExt eq '.bam' or die "[STDERR]: '$bam' does not end in '.bam'\n";
substr($bamDir,0,1) eq '/' or die "[STDERR]: '$bam' is not absolute file name\n";

chdir($bamDir);

my $sortName = $bamName . ".sort";
if ($byReadName) {
	$sortName .= "ByRead";
} else {
	$sortName .= "ByChrom";
}

my $command = "$samtools sort ";
$command .= "-n " if $byReadName;
$command .= " $bam $sortName";

runCommand($command);

my $bamSort = $bamDir . $sortName . $bamExt;

my $BAM2SORT = $bamSort;
if ($replaceOriginal) {
	$BAM2SORT = $bam;
	runCommand("mv $bamSort $bam");
}

runCommand("samtools index $BAM2SORT") unless $byReadName;

print STDERR "done sorting for $bam\n";

sub runCommand {
	my $command = shift;
	
	print STDERR $command, "\n";
	
	!system($command) or die "[STDERR]: can't run $command: $!\n";
}
