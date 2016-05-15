#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use lib '/gpfs/fs121/h/toung/dev';
use GetLineCount;

umask 0007;

my $HOME = '/gpfs/fs121/h/toung';

my $bam;
my $num_threads;
my $tag;
my $gtf_file;
my $frag_bias_correct;
my $cufflinks_path;

$|++;

my $option = GetOptions(
	"bam=s" => \$bam,
	"num_threads=i" => \$num_threads,
	"tag=s" => \$tag,
	"gtf_file=s" => \$gtf_file,
	"frag_bias_correct=s" => \$frag_bias_correct,
	"cufflinks_path=s" => \$cufflinks_path
);

### unsplit options ############################################################
defined $num_threads or die "[STDERR]: define $num_threads\n";
my $options = "-p $num_threads";
-e $gtf_file or die "[STDERR]: $gtf_file doesn't exist\n";
$options .= " -G $gtf_file";
-e $frag_bias_correct or die "[STDERR]: need fasta file $frag_bias_correct\n";
$options .= " --frag-bias-correct $frag_bias_correct";
$options .= " -u"; ### multi-read-correct
################################################################################

### make output directory ######################################################
-e $bam or die "[STDERR]: $bam doesn't exist\n";
my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
$bam_ext eq '.bam' or die "[STDERR]: $bam doesn't end in '.bam'\n";
substr($bam,0,1) eq '/' or die "[STDERR]: $bam not an absolute path\n";
my $exp_dir = $bam_dir . $bam_name . "/cufflinks_$tag";
$exp_dir =~ s/\/aln\//\/exp\// or die "[STDERR]: error creating $exp_dir\n";

runCommand("mkdir -p $exp_dir") unless (-d $exp_dir);

my $command;
$command = $cufflinks_path . "/" if $cufflinks_path;
$command = $command . "cufflinks $options -o $exp_dir $bam";

runCommand($command);

print STDERR "completed_cufflinks:\t$bam\n";
