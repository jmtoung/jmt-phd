#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $fasta;
my $output;
my $blat;
my $index;
my $kmer = 12;
my $stepSize = 5;
my $minIdentity;
my $repMatch = 2253;
my $out = "pslx";

my $options = GetOptions(
	"fasta=s" => \$fasta,
	"output=s" => \$output,
	"blat=s" => \$blat,
	"index=s" => \$index,
	"kmer=i" => \$kmer,
	"stepSize=i" => \$stepSize,
	"minIdentity=i" => \$minIdentity,
	"repMatch=i" => \$repMatch,
	"out=s" => \$out
);

(print STDERR "fasta:\t$fasta\n") && -e $fasta or die "[STDERR]: fasta doesn't exist\n";
(print STDERR "blat:\t$blat\n") && -e $blat or die "[STDERR]: cant find blat\n";
(print STDERR "index:\t$index\n") && -e $index or die "[STDERR]: undefined $index\n";
(print STDERR "kmer:\t$kmer\n") && defined $kmer or die "[STDERR]: undefined kmer\n";
(print STDERR "stepSize:\t$stepSize\n") && defined $stepSize or die "[STDERR]: undefined stepSize\n";
(print STDERR "minIdentity:\t$minIdentity\n") && defined $minIdentity or die "[STDERR]: undefined minIdentity\n";
(print STDERR "repMatch:\t$repMatch\n") && defined $repMatch or die "[STDERR]: undefined repMatch\n";

(print STDERR "output:\t$output\n") && defined $output or die "[STDERR]: output not defined\n";

my $blatCommand = "time $blat -stepSize=$stepSize -minIdentity=$minIdentity -repMatch=$repMatch -noHead -out=$out $index $fasta $output";
print STDERR $blatCommand, "\n";
!system($blatCommand) or die "[STDERR]: can't run $blatCommand: $!\n";

print STDERR "finishedRunning $blatCommand\n";
