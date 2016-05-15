#!/usr/bin/perl -w

use strict;
use Getopt::Long;

### help variable
my $help;
### merged bam file
my $samfile = "/home/jmtoung/Lab/rdd/complete_genomics/ASM_Build36/NA06985/GS06985-1100-36-ASM/GS00392-DNA_D01/ASM/GM06985_pilot_merge/GM06985_pilot_merge_sort.sam";
my $mapfiles = "/home/jmtoung/Lab/rdd/complete_genomics/ASM_Build36/NA06985/GS06985-1100-36-ASM/GS00392-DNA_D01/ASM/GM06985_pilot_1/GM06985_pilot_1.map,/home/jmtoung/Lab/rdd/complete_genomics/ASM_Build36/NA06985/GS06985-1100-36-ASM/GS00392-DNA_D01/ASM/GM06985_pilot_2/GM06985_pilot_2.map";

### open samfile
open(SAM,$samfile);

### this variable keeps track of the read_name of all the alignments we're currently working on
my $prev_readname;
### this hash keeps all the alignments we're working on
my %alignments;

### reference label
my $reference_label = "GM06985_pilot_chr22";

### this is the results/summary hash
my %results;

### add mapfiles to hash
my %mapfiles;
foreach my $mapfile (split(',',$mapfiles)) {
	### get label id for map file
	my @split = split('/',$mapfile);
	my $label_id = pop(@split);
	if ($label_id =~ /(.*).map/) {
		$label_id = $1;
	}
	$mapfiles{$label_id} = $mapfile;
}

### loop through the sam file
while(<SAM>) {
	chomp;
	my @split = split('\t');

	next if $split[1] == 4;
		
	my $readname = $split[0];

	### if prev_readname is undefined, this means we're on the first line, so initalize the varaible and then next;
	if (!defined $prev_readname) {
		### set prev_readname
		$prev_readname = $split[0];
		put_in_hash(\@split,\%alignments);
	### if prev_readname is not undefined, check if its same as previous, if it is, add it to hash
	} elsif ($readname eq $prev_readname) {
		put_in_hash(\@split,\%alignments);
	### if prev_readname is not undefined and if its not same as previous, then process the previous batch and start a new batch
	} else {
		### process previous alignments
		process_batch(\%alignments,\%results, $reference_label, \%mapfiles);
		### clear alignments
		%alignments = ();
		### set prev_readname to current readname
		$prev_readname = $split[0];
		put_in_hash(\@split,\%alignments);
	}
}


process_batch(\%alignments, \%results, $reference_label, \%mapfiles);

$DB::single = 1;
print "hi\n";



sub process_batch {
	my $alignments = shift(@_);
	my $results = shift(@_);
	my $reference_label = shift(@_);
	my $mapfiles = shift(@_);
	
	### this is a temporary results hash to tally things up before you put in final results hash
	my %temp_results;
	
	### this variable tells you what group the read belongs in: reference or individual
	my $group;
	### convert coordinates
	foreach my $number (keys %{$alignments{"alignments"}}) {

		### if the group field (label) is the reference label, then the group is 'ref', else it's 'ind'
		my $label = $$alignments{"alignments"}{$number}{"RG"};
		if ($label eq $reference_label) {
			$group = "ref";
		} else {
			$group = "ind";
		}
		

		### this means there's no alignment
		if ($$alignments{"alignments"}{$number}{"flag"} == 4) {
			next;
		} else {
			my $chrom = $$alignments{"alignments"}{$number}{"chrom"};
			my $position = $$alignments{"alignments"}{$number}{"position"};
			
			### convert coordinates
			if ($group eq "ind") {
				
				my $mapfile = $$mapfiles{$label};
				$position = convert_position($position,$mapfile);
			}
			my $num_mismatch = num_mismatch($$alignments{"alignments"}{$number}{"MD"});
			
			$temp_results{$chrom}{$position}{$group}{$num_mismatch}++;
		}
	}
	
	### tally up results for the results hash
	foreach my $chrom (keys %temp_results) {
		### number of groups withat have this alignment
		foreach my $position (keys %{$temp_results{$chrom}}) {
			my @groups = sort keys %{$temp_results{$chrom}{$position}};
			### if there are both ref and ind alignments
			if ($#groups == 1) {
				### difference in number of mismatches
				my $diff = scalar(keys %{$temp_results{$chrom}{$position}{"ref"}}) - scalar(keys %{$temp_results{$chrom}{$position}{"ind"}});
				$$results{"venn"}{"both"}{$diff}++;
			} elsif ($groups[0] eq 'ref') {
				my $num_mismatch = keys %{$temp_results{$chrom}{$position}{"ref"}};
				$$results{"venn"}{"ref"}{$num_mismatch}++;
			} elsif ($groups[0] eq 'ind') {
				my $num_mismatch = keys %{$temp_results{$chrom}{$position}{"ind"}};
				$$results{"venn"}{"ind"}{$num_mismatch}++;
			}
						
		}
		
		
	}
}

sub convert_position {
	my $position = shift(@_);
	my $mapfile = shift(@_);
	
        my $awk = "awk '" . '$4 <=' . $position . "' $mapfile | tail --lines=1 |";
        open(AWK,$awk) or die "couldn't fork: $!\n";

	my @map_values;        
        while(<AWK>) {
		@map_values = split('\t');
	}
	
	### if it's an insertion 
	if ($map_values[2] < $map_values[1]) {
		return -1;
	} else {
		return $position - $map_values[3] + $map_values[1];
	}
}

### this subroutine gives you the number of mismatches (essentially the number of A/C/T/G)
sub num_mismatch {
	my $mismatch_string = shift(@_);
	
	### split the string
	my @split = split('',$mismatch_string);
	
	### num of mismatch
	my $num_mismatch = 0;
	foreach my $letter (@split) {
		if ($letter eq 'A' || $letter eq 'C' || $letter eq 'G' || $letter eq 'T') {
			$num_mismatch++;
		}
	}
	return $num_mismatch;
}

sub put_in_hash {
	my $split = shift(@_);
	my $alignments = shift(@_);
	my $labels = shift(@_);
	
	my $flag = $split->[1];
	my $chrom = $split->[2];
	my $position = $split->[3];
	my @optional_fields = @$split[11..(scalar(@$split) - 1)];
	
	
	### the number of the alignment.
	my $number;
	### check and see if hash has already been initalized
	unless (defined $$alignments{"number"}) {
		### if it hasn't been defined, set the number to 1. the number keeps track of all the alignments 
		$number = $$alignments{"number"} = 1;
	} else {
		$number = ++$$alignments{"number"};
	}
	
	### now put all the info for this alignment into the hash
	$$alignments{"alignments"}{$number}{"flag"} = $flag;
	$$alignments{"alignments"}{$number}{"chrom"} = $chrom;
	$$alignments{"alignments"}{$number}{"position"} = $position;
	
	### parse the optional fields and put it into hash
	foreach my $optional_field (@optional_fields) {
		my @split = split(':',$optional_field);
		$$alignments{"alignments"}{$number}{$split[0]} = $split[$#split];
	}
	
	

}


