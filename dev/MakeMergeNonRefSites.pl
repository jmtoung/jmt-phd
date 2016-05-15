#!/usr/bin/perl -w

use strict;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Database;
use File::Basename;
use ComplementBase;
use Getopt::Long;

my $HOME = "/home/jmtoung/Lab";

my $collection;
my $alnDB = "$HOME/aln/alnDB";
my $tag = ".Qual20NonRefExCompPoly.txt";
my $min_total_reads = 10;
my $min_rdd_level = 0;
my $single_nuc_changes = 1;

my $options = GetOptions(
	"collection=s" => \$collection,
	"alnDB=s" => \$alnDB,
	"tag=s" => \$tag,
	"min_total_reads=i" => \$min_total_reads,
	"min_rdd_level=s" => \$min_rdd_level,
	"single_nuc_changes=i" => \$single_nuc_changes
);

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
print STDERR "option_min_total_reads:", "\t", $min_total_reads, "\n";
print STDERR "option_min_rdd_level:", "\t", $min_rdd_level, "\n" if $min_rdd_level;
print STDERR "option_single_nuc_changes:", "\t", $single_nuc_changes, "\n\n";

### REPLACE DOTS (.) IN $tag WITH (\.) #########################################
defined $tag or die "[STDERR]: tag not defined\n";
$tag =~ s/\./\\./g;
################################################################################

### LOAD THE PILEUP DIRECTORIES FOR THE bamIDs #################################
my @pileup_dir;
foreach my $bamID (@bamID) {
	my $bam = Database->new($alnDB)->lookup(0,$bamID,1);
	-e $bam or die "[STDERR]: $bam doesn't exist: $!\n";
	my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
	my $pileup_dir = $bam_dir . $bam_name . ".unique.pileup";
	-d $pileup_dir or die "[STDERR]: $pileup_dir doesn't exist: $!\n";
	push(@pileup_dir,$pileup_dir);
}
################################################################################

### LOAD NONREF_SITES FILES FOR EACH PILEUP DIRECTORY BY CHROMOSOME ############
my %nonref_sites;
foreach my $pileup_dir (@pileup_dir) {
	my $find_statement = "find $pileup_dir -name \\*$tag |";
	print STDERR "find_statement", "\t", $find_statement, "\n";

	open(FIND,$find_statement) or die "[STDERR]: can't fork $find_statement: $!\n";
	while(<FIND>) {
		chomp;
		print STDERR "...found...", "\t", $_, "\n";
		my ($nonref_name, $nonref_dir, $nonref_ext) = fileparse($_,$tag);
		my @nonref_name = split('\.',$nonref_name);
		my $chrom = $nonref_name[-1];
		$nonref_sites{$chrom}{$_}++;
	}
}
################################################################################

### LOOP THROUGH EACH CHROMOSOME ###############################################
my @FILES;
foreach my $chrom (sort keys %nonref_sites) {
	my @files = sort keys %{$nonref_sites{$chrom}};
	print STDERR "chromosome:", "\t", $chrom, "\t", "num_files:", "\t", scalar(@files), "\n";

	my $files = join(",",@files);
	push(@FILES,$files);
}

my $FILES = join("\n",@FILES);
my $num_files = scalar(@FILES);

print STDOUT <<"END";
#\$ -cwd ### use current directory
#\$ -S /usr/bin/perl ### program to execute script
#\$ -M toung\@mail.med.upenn.edu ### email address
##\$ -m ea ### mail is to be sent at abort and end time
#\$ -t 1-$num_files ### array job #
#\$ -V ### use current environment variables
##\$ -pe DJ 2 
##\$ -l mem_free=4G

use strict;

my \$ID = \$ENV{SGE_TASK_ID} - 1;

my \@files = qw(
$FILES
);

print STDERR "going_to_do\\t\$files[\$ID]\\n";

system("perl $HOME/dev/MergeNonRefSites.pl --collection $collection --files \$files[\$ID] --min_total_reads $min_total_reads --min_rdd_level $min_rdd_level --single_nuc_changes $single_nuc_changes");

print STDERR "completed\\t\$files[\$ID]\\n";

END

