#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use File::Basename;
use Bio::SeqIO;
use Bio::Seq;

my $fastq;
my $tag;

my $options = GetOptions(
	"fastq=s" => \$fastq,
	"tag=s" => \$tag,
);

my($fastq_name,$fastq_dir,$fastq_ext) = fileparse($fastq,'\.fastq(\.gz)?');

my $command;
if ($fastq_ext =~ /\.gz/) { $command = "gunzip -c $fastq |"; }
elsif ($fastq_ext eq '.fastq') { $command = $fastq; }
else { die "[STDERR]: can't open $fastq\n"; }


$fastq_ext =~ s/\.gz//;
my $output = $fastq_dir . $fastq_name . "_" . $tag . $fastq_ext;
open(OUTPUT,">".$output) or die "[STDERR]: can't open $output: $!\n";

open(FASTQ,$command) or die "[STDERR]: can't fork $command: $!\n";

while(my $line1 = <FASTQ>) {
        defined (my $line2 = <FASTQ>) or die "[STDERR]: line2 (sequence) undefined\n";
        defined (my $line3 = <FASTQ>) or die "[STDERR]: line3 (read_name) undefined\n";
        defined (my $line4 = <FASTQ>) or die "[STDERR]: line4 (qual scores) undefined\n";
        chomp ($line1, $line2, $line3, $line4);

	my $dna_obj = Bio::Seq->new(-seq => $line2, -display_id => "", -alphabet => "dna");
	$line2 = $dna_obj->revcom->seq;

	print OUTPUT $line1, "\n", $line2, "\n", $line3, "\n", scalar(reverse($line4)), "\n";
}
