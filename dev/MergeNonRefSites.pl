#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Database;
use File::Basename;
use ComplementBase;

my $HOME = "/home/jmtoung/Lab";

my $collection;
my $files;
my $min_total_reads = 10;
my $min_rdd_level = 0;
my $single_nuc_changes = 1;

my $options = GetOptions(
	"collection=s" => \$collection,
	"files=s" => \$files,
	"min_total_reads=i" => \$min_total_reads,
	"min_rdd_level=s" => \$min_rdd_level,
	"single_nuc_changes=i" => \$single_nuc_changes
);

my @files = split(',',$files);

$files = join(" ",@files);

### MAKE CAT STATEMENT #########################################################
my $cat = "cat $files | perl $HOME/dev/ConvertToForwardStrand.pl |";
if ($min_total_reads) { $cat .= " awk '(\$5>=$min_total_reads)' |"; }
if ($min_rdd_level) { $cat .= " awk '(\$7/\$5>=$min_rdd_level)' |"; }
if ($single_nuc_changes) { $cat .= " awk '(\$4~/^[ACGT]\$/ && \$6~/^[ACGT]\$/\)' |"; }
else { $cat .= " awk '(\$4~/^[ACGT]\$/ && \$6~/[ACGT][ACGT]\.\*/\)' |"; }
### sort by chrom, position, ref base, nonref base, strand, alnID 
### DONT sort by strand first because some might be strand-specific
$cat .= " sort +0 -1 +1n -2 +3 -4 +5 -6 +2 -3 +8n -9 |";
################################################################################

print STDERR "cat command is:\t$cat\n";

open(CAT,$cat) or die "[STDERR]: can't fork $cat: $!\n";

my %current;
while(<CAT>) {
	chomp;
	my @split = split('\t');

	if (is_equal(\@split,\%current)) {
		add_to_current(\@split,\%current);
	} else {
		print_current(\%current);
		%current = ();
		add_to_current(\@split,\%current);
	}
}
close(CAT);

print_current(\%current);

##################################################################################

sub print_current {
	my $current = shift;

	print $current->{'chrom'}, "\t";
	print $current->{'position'}, "\t";
	print $current->{'strands'}, "\t";
	print $current->{'ref_base'}, "\t";
	print $current->{'ref_base_count'}, "\t";
	print $current->{'nonref_base'}, "\t";
	print $current->{'nonref_base_count'}, "\t";
	print $current->{'poly'}, "\t";

	my @bamID = sort {$a <=> $b} keys %{$current->{'bamID'}};
	my $bamID = join(",",@bamID);
	print "{", $bamID, "}", "\t";
	print scalar(@bamID), "\t";
	print $collection, "\n";

}

sub is_equal {
	my($split,$current) = @_;

	### if hash is empty (like for the first line), return equal
	return 1 if !exists $current->{'chrom'};

	return 0 unless $current->{'chrom'} eq $split->[0];
	return 0 unless $current->{'position'} == $split->[1];
	return 0 unless $current->{'strands'} eq $split->[2];
	return 0 unless $current->{'ref_base'} eq $split->[3];
	return 0 unless $current->{'nonref_base'} eq $split->[5];

	return 1;
}

sub add_to_current {
	my($split,$current) = @_;

	$current->{'chrom'} = $split->[0];
	$current->{'position'} = $split->[1];
	$current->{'strands'} = $split->[2];
	$current->{'ref_base'} = $split->[3];
	$current->{'nonref_base'} = $split->[5];

	$current->{'ref_base_count'} += $split->[4];
	$current->{'nonref_base_count'} += $split->[6];

	if (!exists $current->{'poly'}) { $current->{'poly'} = $split->[7]; } 
	else { $current->{'poly'} eq $split->[7] or die "[STDERR]: mismatch in poly ", join("\t",@{$split}), "\n"; }

	$current->{'bamID'}{$split->[8]}++;
}
