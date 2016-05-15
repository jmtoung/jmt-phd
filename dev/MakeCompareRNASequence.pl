#!/usr/bin/perl -w

use strict;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use List::Util qw(max min);
use Database;
use File::Basename;
use ComplementBase;
use Getopt::Long;

### change these if necessary
my $HOME = "/ifs/h/toung";
my $chrom_sizes = "/ifs/h/toung/database/chrom_sizes_gsnap.GM12750_ADAR_KD_trimlowqual.def.unique.bam.txt";
###

my $collection;
my $alnDB = "/ifs/h/toung/aln/alnDB";
my $tag = ".pileup";
my $minqual = 20;
my $min_total_reads = 10;
my $max_size = 10000000;

my $options = GetOptions(
	"collection=s" => \$collection,
	"alnDB=s" => \$alnDB,
	"tag=s" => \$tag,
	"minqual=i" => \$minqual,
	"min_total_reads=i" => \$min_total_reads,
	"max_size=i" => \$max_size
);

### LOAD CHROM SIZES ###########################################################
open(CHROM_SIZES,$chrom_sizes) or die "[STDERR]: can't open $chrom_sizes: $!\n";
my %chrom_sizes;
while(<CHROM_SIZES>) {
	chomp;
	my @split = split('\t');
	$chrom_sizes{$split[0]} = $split[1];
}

### LOAD THE bamIDs FROM THE COLLECTIONS FILE ##################################
my @bamID;
my $collection_file = "/ifs/h/toung/database/collectionsDB/" . $collection;
my ($collection_name, $collection_dir, $collection_ext) = fileparse($collection_file,'');
open(COLLECTION,$collection_file) or die "[STDERR]: can't open $collection_file: $!\n";
while(<COLLECTION>) {
	chomp;
	my @split = split('\t');
	my $bamID = $split[1];
	push(@bamID,$bamID);
}
close COLLECTION;
################################################################################

print STDERR "option_collection:", "\t", $collection, "\n";
print STDERR "option_alnDB:", "\t", $alnDB, "\n";
print STDERR "option_tag:", "\t", $tag, "\n";
print STDERR "option_minqual:", "\t", $minqual, "\n";
print STDERR "option_min_total_reads:", "\t", $min_total_reads, "\n\n";

### REPLACE DOTS (.) IN $tag WITH (\.) #########################################
defined $tag or die "[STDERR]: tag not defined\n";
$tag =~ s/\./\\./g;
################################################################################

### LOAD THE PILEUP DIRECTORIES FOR THE bamIDs #################################
my %pileup_dir;
foreach my $bamID (@bamID) {
	my $bam = Database->new($alnDB)->lookup(0,$bamID,1);
	-e $bam or die "[STDERR]: $bam doesn't exist: $!\n";
	my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
	my $pileup_dir = $bam_dir . $bam_name . ".unique.pileup";
	-d $pileup_dir or die "[STDERR]: $pileup_dir doesn't exist: $!\n";
	$pileup_dir{$bamID} = $pileup_dir;
}
################################################################################

### LOAD PILEUP FILES FOR EACH PILEUP DIRECTORY BY CHROMOSOME ##################
my %pileup_files;
foreach my $bamID (keys %pileup_dir) {
	my $find_statement = "find $pileup_dir{$bamID} -name \\*$tag -type f | grep '.*\\.[0-9]\\{1,2\\}\\.pileup' |";
	print STDERR "find_statement", "\t", $find_statement, "\n";

	open(FIND,$find_statement) or die "[STDERR]: can't fork $find_statement: $!\n";
	while(<FIND>) {
		chomp;
		print STDERR "...found...", "\t", $_, "\n";
		my ($pileup_name, $pileup_dir, $pileup_ext) = fileparse($_,$tag);
		my @pileup_name = split('\.',$pileup_name);
		my $chrom = $pileup_name[-1];
		$_ = $bamID . ":" . $_;
		$pileup_files{$chrom}{$_}++;
	}
}
################################################################################

### LOOP THROUGH EACH CHROMOSOME ###############################################
my @FILES;
my @START;
my @END;

foreach my $chrom (sort keys %pileup_files) {

	my @files = sort keys %{$pileup_files{$chrom}};
	my $files = join(",",@files);

	print STDERR "chromosome:", "\t", $chrom, "\t", "num_files:", "\t", scalar(@files), "\n";

        ### divy up each chromosome according to max_size & length of chrom ####
        my $START = 1;
        my $END = min($chrom_sizes{$chrom}, $max_size);
        do {
		push(@START,$START);
		push(@END,$END);
		push(@FILES,$files);

                $START = ++$END;
                $END = min($chrom_sizes{$chrom}, $START + $max_size - 1);
        } until ($START > $chrom_sizes{$chrom});
        ########################################################################
}

my $FILES = join("\n",@FILES);
my $START = join("\n",@START);
my $END = join("\n",@END);

my $num_files = scalar(@FILES);

print STDOUT <<"END";
#\$ -cwd ### use current directory
#\$ -S /usr/bin/perl ### program to execute script
#\$ -M toung\@mail.med.upenn.edu ### email address
##\$ -m ea ### mail is to be sent at abort and end time
#\$ -t 1-$num_files ### array job #
#\$ -V ### use current environment variables
#\$ -pe DJ 2
#\$ -l mem_free=2G

use strict;

my \$ID = \$ENV{SGE_TASK_ID} - 1;

my \@files = qw(
$FILES
);

my \@start = qw(
$START
);

my \@end = qw(
$END
);

print STDERR "perl $HOME/dev/CompareRNASequence.pl --pileup \$files[\$ID] --minqual $minqual --min_total_reads $min_total_reads --start_position \$start[\$ID] --end_position \$end[\$ID]\\n";

system("perl $HOME/dev/CompareRNASequence.pl --pileup \$files[\$ID] --minqual $minqual --min_total_reads $min_total_reads --start_position \$start[\$ID] --end_position \$end[\$ID]");

print STDERR "completed\\t\$files[\$ID]\\n";

END

