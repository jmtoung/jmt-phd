#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Database;
use File::Basename;
use ComplementBase;
use Tie::IxHash;

my $HOME = "/ifs/h/toung";

my $file; 
my $min_rdd_reads;
my $min_individuals;

my $options = GetOptions(
	"file=s" => \$file,
	"min_rdd_reads=s" => \$min_rdd_reads,
	"min_individuals=s" => \$min_individuals
);

$|++;

-e $file or die "[STDERR]: file $file doesn't exist:\n";

print STDOUT "file:\t$file\n";
print STDOUT "min_rdd_reads:\t$min_rdd_reads\n" if defined $min_rdd_reads;
print STDOUT "min_individuals\t$min_individuals\n" if defined $min_individuals;

### only one of the two (min_total_reads & min_individuals) can be defined #####
if ($min_rdd_reads || $min_individuals) { die "[STDERR]: one and only one of min_rdd_reads and min_individuals need to be defined\n" if $min_rdd_reads && $min_individuals;
} else { die "[STDERR]: one and only one of min_rdd_reads and min_individuals need to be defined\n"; }
################################################################################

### DEFINED THRESHOLDS FOR THE SUBLISTS ########################################
my %filters;
if (defined $min_individuals) {
	if ($min_individuals =~ /-/) {
		my ($start,$end) = split('-',$min_individuals);
		$filters{'min_individuals'} = [$start..$end];
	} elsif ($min_individuals =~ /,/) {
		$filters{'min_individuals'} = [split(',',$min_individuals)] if $min_individuals;
	}
} elsif (defined $min_rdd_reads) {
	if ($min_rdd_reads =~ /-/) {
		my ($start,$end) = split('-',$min_rdd_reads);
		$filters{'min_rdd_reads'} = [$start..$end];
	} elsif ($min_rdd_reads =~ /,/) {
		$filters{'min_rdd_reads'} = [split(',',$min_rdd_reads)] if $min_rdd_reads;
	}
}

### PARSE FILE NAME FOR OUTPUTS ################################################
my ($file_name,$file_dir,$file_ext) = fileparse($file,'\.txt');

### MAKE SUBLISTS ##############################################################
my %sublists;
tie %sublists, 'Tie::IxHash'; ### ordered hash
$sublists{'total_count'} = "cat $file"; ### no filter
$sublists{'total_count_rm_poly'} = "awk '(\$8~/{-}/)' $file"; ### remove poly
### min_total_reads OR min_individuals
foreach my $filter_name (keys %filters) {
	foreach my $threshold (@{$filters{$filter_name}}) {
		my $sublist_name = "total_count_rm_poly_" . $filter_name . "_" . $threshold;
		$sublists{$sublist_name} = "awk '(\$8~/{-}/ &&";
		if ($filter_name eq 'min_rdd_reads') { $sublists{$sublist_name} .= " \$7>=$threshold)'"; } 
		elsif ($filter_name eq 'min_individuals') { $sublists{$sublist_name} .= " \$10>=$threshold)'"; }
		$sublists{$sublist_name} .= " $file";
	}
}

### DETERMINE COLLECTION NAME AND BAM IDS IN COLLECTION ########################
my $collection = get_collection($sublists{'total_count'});
my $collection_file = $HOME . "/database/collectionsDB/" . $collection;
open(COLLECTIONS,$collection_file) or die "[STDERR]: can't open collection file $collection_file: $!\n";
my @collections;
while(<COLLECTIONS>) {
	chomp;
	my @split = split('\t');
	push(@collections,$split[1]);
}
@collections = sort {$a <=> $b} @collections;

### CALCULATE COUNTS  ###########################################################
my $counts_output = $file_dir . $file_name . ".counts";
open(COUNTS,">".$counts_output) or die "[STDERR]: can't open $counts_output: $!\n";
foreach my $sublist_name (keys %sublists) {
	my $command = $sublists{$sublist_name} . " | awk '{print \$4,\$6}' | sort | uniq -c | awk 'BEGIN {OFS=\"\\t\"} {print \$1,\$2,\$3}' |";
	print STDOUT "count_command:\t$sublist_name\t$command\n";

	open(COMMAND,$command) or die "[STDERR]: can't fork $command: $!\n";
	while(<COMMAND>) {
		chomp;
		my @split = split('\t');
		print COUNTS $collection, "\t", $sublist_name, "\t", $split[1], "\t", $split[2], "\t", $split[0], "\n";
	} 
	close(COMMAND);
}
close(COUNTS);

### CALCULATE SET DISTRIBUTION #################################################
my $set_output = $file_dir . $file_name . ".set";
open(SET_OUTPUT,">".$set_output) or die "[STDERR]: can't open $set_output: $!\n";
my $binary_set_output = $file_dir . $file_name . ".bset";
open(BINARY_SET_OUTPUT,">".$binary_set_output) or die "[STDERR]: can't open $binary_set_output: $!\n";

foreach my $sublist_name (keys %sublists) {
	my $command = $sublists{$sublist_name} . " | awk '{print \$4,\$6,\$9}' | sort +2 -3 +0 -1 +1 -2 | uniq -c | awk 'BEGIN {OFS=\"\\t\"} {print \$1,\$2,\$3,\$4}' |";
	print STDOUT "set_command:\t$sublist_name\t$command\n";

	open(COMMAND,$command) or die "[STDERR]: can't fork $command: $!\n";

	my %binary_sets; ### stores total counts for each pairwise in the group
	while(<COMMAND>) {
		chomp;
		my @split = split('\t');

		### print to the total set output list
		print SET_OUTPUT $collection, "\t", $sublist_name, "\t", $split[0], "\t", $split[1], "\t", $split[2], "\t", $split[3], "\n";

		### parse the list of bamIDs for binary set computation
		(my $bamIDs = $split[3]) =~ s/[{}]//g;
		my @bamIDs = split(',',$bamIDs);
		@bamIDs = sort {$a<=>$b} @bamIDs;
		
		### store counts for each binary combination of the set in @bamIDs
		foreach my $bamID_1 (@bamIDs) {
			foreach my $bamID_2 (@bamIDs) {
				next if $bamID_1 == $bamID_2;
				$binary_sets{$bamID_1}{$bamID_2}{$split[1]}{$split[2]} += $split[0];
			}
		} 

	}
	close(COMMAND);

	### print out counts for each binary combination
	foreach my $bamID_1 (sort {$a<=>$b} keys %binary_sets) {
		foreach my $bamID_2 (sort {$a<=>$b} keys %{$binary_sets{$bamID_1}}) {		
			foreach my $ref_base (sort keys %{$binary_sets{$bamID_1}{$bamID_2}}) {
				foreach my $rdd_base (sort keys %{$binary_sets{$bamID_1}{$bamID_2}{$ref_base}}) {
					print BINARY_SET_OUTPUT $collection, "\t", $sublist_name, "\t";
					print BINARY_SET_OUTPUT $binary_sets{$bamID_1}{$bamID_2}{$ref_base}{$rdd_base}, "\t";
					print BINARY_SET_OUTPUT $ref_base, "\t", $rdd_base, "\t";
					print BINARY_SET_OUTPUT "{", $bamID_1, ",", $bamID_2, "}\n";

				}
			}
		}
	}
}
close(SET_OUTPUT);
close(BINARY_SET_OUTPUT); 

################################################################################
sub get_collection {
	my $command = shift;

	$command .= " | head -1 | awk '{print \$11}'";
	
	my $collection = `$command`;
	chomp $collection;

	return $collection;
}
