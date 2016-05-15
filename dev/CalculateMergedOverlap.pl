#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Database;
use File::Basename;
use ComplementBase;
use Tie::IxHash;
use List::Util qw(max min);

my $HOME = "/ifs/h/toung";

my $files;
my $tag;
my $filters;

my $options = GetOptions(
	"files=s" => \$files,
	"tag=s" => \$tag,
	"filters=s" => \$filters
);

my @files = split(',',$files);
foreach my $file (@files) { print "file:\t$file\n"; }
my @filters = split('Z',$filters);
foreach my $filter (@filters) { print "filter:\t$filter\n"; }
print "tag:\t$tag\n";
scalar(@files) == scalar(@filters) or die "[STDERR]: each file must get a filter\n";
$files = join(" ",@files);

my %FILTERS;
tie %FILTERS, 'Tie::IxHash';
$FILTERS{'total'} = "cat $files";
$FILTERS{'total_rm_poly'} = "awk '(\$8~/{-}/)' $files";

### add other complicated filters
my %FILTERS_COMP;
tie %FILTERS_COMP, 'Tie::IxHash';

foreach my $file (@files) {
	### determine collection for file
	my $collection = get_collection($file);

	### determine the type of filter
	my $filter = shift(@filters);
	my ($filter_type,$filter_thresholds) = split(':',$filter);

	### determine column of filter
	my $filter_column;
	if ($filter_type eq 'min_rdd_reads') { $filter_column = 7; }
	elsif ($filter_type eq 'min_individuals') { $filter_column = 10; }
	else { die "[STDERR]: filter type must be 'min_rdd_reads' or 'min_individuals'\n"; }

	### determine filter thresholds
	my @filter_thresholds;
	if ($filter_thresholds =~ /([0-9]+)-([0-9]+)/) { @filter_thresholds = ($1..$2); }
	elsif ($filter_thresholds =~ /,/) {
		@filter_thresholds = split(',',$filter_thresholds);
	} elsif ($filter_thresholds =~ /[0-9]+/) {
		@filter_thresholds = ($filter_thresholds);
	} else { die "[STDERR]: invalid $filter_thresholds\n"; }

	### if there are currently no filters, add them to %FILTERS_COMP (prefaced by awk '($8~/{-}/)'
	if (scalar(keys %FILTERS_COMP) == 0) {
		foreach my $filter_threshold (@filter_thresholds) {
			my $filter_name = $collection . "_" . $filter_type . "_" . $filter_threshold;
			$FILTERS_COMP{$filter_name} = "awk '(((\$11~/$collection/ && \$$filter_column >= $filter_threshold)";
		}
	} else {
	### add new filters to the end of each pre-existing
		my %FILTERS_COMP_ADDL;
		tie %FILTERS_COMP_ADDL, 'Tie::IxHash';
		foreach my $filter_name (keys %FILTERS_COMP) {
			foreach my $filter_threshold (@filter_thresholds) {
				my $filter_name_addl = $filter_name . "_" . $collection . "_" . $filter_type . "_" . $filter_threshold;
				$FILTERS_COMP_ADDL{$filter_name_addl} = $FILTERS_COMP{$filter_name} . " || (\$11~/$collection/ && \$$filter_column >= $filter_threshold)";
			}
		}
		%FILTERS_COMP = ();
		%FILTERS_COMP = %FILTERS_COMP_ADDL;
	}
}

### complete the awk statements
foreach my $filter_name (keys %FILTERS_COMP) {
	$FILTERS{$filter_name} = $FILTERS_COMP{$filter_name} . "))' $files";
}

my $output_counts = "MergeMergedNonRefSites_$tag.counts";
print "output_counts:\t$output_counts\n";
my $output_set = "MergeMergedNonRefSites_$tag.set";
print "output_set:\t$output_set\n";
open(OUTPUT_COUNTS,">".$output_counts) or die "[STDERR]: can't open $output_counts: $!\n";
open(OUTPUT_SET,">".$output_set) or die "[STDERR]: can't open $output_set: $!\n";

foreach my $filter_name (keys %FILTERS) {
	### sort by chrom, position, ref_base, rdd_base, strand, collection
	my $command = $FILTERS{$filter_name} . " | sort +0 -1 +1n -2 +3 -4 +5 -6 +2 -3 +10 -11 |";

	print STDOUT "\ncommand\t$filter_name\t", $command, "\n";
	open(COMMAND,$command) or die "[STDERR]: can't fork $command: $!\n";

	my $output = "MergeMergedNonRefSites_$tag-$filter_name.txt";
	print STDOUT "output\t$filter_name\t$output\n";
	open(OUTPUT,">".$output) or die "[STDERR]: can't open $output: $!\n";
	select(OUTPUT);

	my (%counts,%coll_group); ### tally up data
	my %current;
	while(<COMMAND>) {
		chomp;
		my @split = split('\t');
	
		if (is_equal(\@split,\%current)) {
			add_to_current(\@split,\%current);
		} else {
			print_current(\%current,\%counts,\%coll_group);
			%current = ();
			add_to_current(\@split,\%current);
		}
	}
	close(COMMAND);

	print_current(\%current,\%counts,\%coll_group);
	close(OUTPUT);
	
	foreach my $ref_base (keys %counts) {
		foreach my $rdd_base (keys %{$counts{$ref_base}}) {
			print OUTPUT_COUNTS "$filter_name", "\t", $ref_base, "\t", $rdd_base, "\t", $counts{$ref_base}{$rdd_base}, "\n";
		}
	}

	foreach my $collection_group (keys %coll_group) {
		print OUTPUT_SET $filter_name, "\t", $collection_group, "\t", $coll_group{$collection_group}, "\n";
	}
}

##################################################################################

sub print_current {
	my $current = shift;
	my $counts = shift;
	my $coll_group = shift;

	print $current->{'chrom'}, "\t";
	print $current->{'position'}, "\t";
	print $current->{'strands'}, "\t";
	print $current->{'ref_base'}, "\t";
	print $current->{'ref_base_count'}, "\t";
	print $current->{'nonref_base'}, "\t";
	print $current->{'nonref_base_count'}, "\t";
	print $current->{'poly'}, "\t";

	print "{", join(",",@{$current->{'bamIDset'}}), "}", "\t";
	print min(@{$current->{'num_bamIDset'}}), "\t";
	my $collection_group = join(",",@{$current->{'collection'}});
	print "{", $collection_group, "}", "\t";
	print scalar(@{$current->{'collection'}}), "\n";
	
	$counts->{$current->{'ref_base'}}{$current->{'nonref_base'}}++;
	$coll_group->{$collection_group}++;
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

	push(@{$current->{'bamIDset'}},$split->[8]);
	push(@{$current->{'num_bamIDset'}},$split->[9]);
	push(@{$current->{'collection'}},$split->[10]);
}

sub get_collection {
        my $file = shift;

	-e $file or die "[STDERR]: $file doesn't exist\n";

	my $command = "cat $file | head -1 | awk '{print \$11}'";

        my $collection = `$command`;
        chomp $collection;

        return $collection;
}
