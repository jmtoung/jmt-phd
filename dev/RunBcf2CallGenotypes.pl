#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $bam;
my $index;
my $tag;
my $options;

my $opt = GetOptions(
	"bam=s" => \$bam,
	"index=s" => \$index,
	"tag=s" => \$tag,
	"options=s" => \$options
);

my @options_parsed;
if ($options) {
	my @options = split(';',$options);
	foreach my $options (@options) {
		my @options_split = split(':',$options);
		push(@options_parsed,@options_split);
	}
}
my $options_parsed = join(" ",@options_parsed);

-e $bam or die "[STDERR]: bam $bam doesn't exist: $!\n";
print STDERR "bam:\t$bam\n";

-e $index or die "[STDERR]: index $index doesn't exist: $!\n";
print STDERR "index:\t$index\n";

defined $tag or die "[STDERR]: please define $tag\n";
print STDERR "tag:\t$tag\n";

print "options:\t$options_parsed\n";

my $output = $bam;
$output =~ s/\.bam/\.$tag\.bcf/;

my $bcf = "samtools mpileup";
$bcf .= " " . $options_parsed if $options;
$bcf .= " -ugf $index $bam | bcftools view -bNcg - > $output";
print STDERR "bcf command is:\n";
print STDERR $bcf, "\n";

!system($bcf) or die "[STDERR]: can't run $bcf: $!\n";

print STDERR "completed\t$bam\n";
