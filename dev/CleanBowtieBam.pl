#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

### bam file
my $bam;
### maximum edit distance (mismatches in bowtie since bowtie doesn't allow indels) allowed
my $max_edit_distance;
### maximum strata (specification only in bowtie)
my $max_strata;

my $result = GetOptions("bam=s" => \$bam, "max_edit_distance=s" => \$max_edit_distance, "max_strata=s" => \$max_strata);

################################################################################
### This script will do four things to a Bowtie alignment file.
### (1) It will filter all alignments that have an edit distance (number of mismatches) greater than $max_edit_distance OR have a strata greater than $max_strata.
### The parameters for $max_edit_distance and $max_strata are fully customizable.
### (2) The script will collapse  all alignments that have the same exact info except they align to different reference sequences. This will occur if the target sequence
### is represented by more than one SS index (in addition to human genome sequence). This filtering step will necessarily occur.
### (3) For each read, the script will place info of other alignments to the end of the string in an 'OT' tag (a la BWA alignment's 'XA' tag). The 'OT' tag is (chr,(+/-)pos,CIGAR,NM;)*.
### (4) This script will change the map qual column such that it corresponds to total number of alignments this read has!

### check that $bam ends in '.bam' and is an absolute file name ################
my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
$bam_ext eq '.bam' or die "[STDERR]: '$bam' does not end in '.bam'\n";
substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$bam' is not absolute file name\n";
chdir($bam_dir) or die "[STDERR]: cannot chdir '$bam_dir': $!\n";
################################################################################

### open output $CLEAN_SAM #####################################################
my $clean_sam = $bam_dir . $bam_name . "_clean.sam";
open(my $CLEAN_SAM,">".$clean_sam) or die "[STDERR]: cannot open $clean_sam: $!\n";
################################################################################

### open the bam file to read using `samtools view` ############################
my $BAM = "samtools view -h " . $bam . " |";
open(BAM,$BAM) or die "[STDERR]: couldn't fork: $!\n";
################################################################################

### loop through BAM, group alignments w/same read name and filter alignments ##
my $counter = 0;
my %aln; ### hash that contains all the alignments
my $read_name = ""; ### keep track of the current read_name
print ">>> Filtering Reads <<<\n";
while(<BAM>) {

	### if it's the header line, just print and leave ######################
	if (/^@/) { 
		print {$CLEAN_SAM} $_; 
		next; 
	}
	
	chomp;
	my @split = split('\t');

	### if we're on the first alignment, then $read_namne is '', so initialize then move on
	if ($read_name eq '') { 
		add_alignment(\@split,\%aln);
		$read_name = $split[0];
		next;
	}

	########################################################################		
	### if $split[0] (current read_name) is equal to the previous read_name ($read_name), then...add alignment!
	### if not, we'll all previous process alignments and then clear variables and add current alignment 
	if ($split[0] eq $read_name) {
		add_alignment(\@split,\%aln);
	} else {
		### collapse alignments that align to the identical positions, identical flags, everything, just align to different ss (or to genomic)	
		collapse_alignments_same_ss(\%aln);
		
		### filter out bad alignments and print
		filter_and_print_alignments(\%aln);

		print "... Finished $counter Reads ...\n" if (++$counter % 1000000) == 0;

		### clear alignments
		%aln = ();
		### add new read to alignments
		add_alignment(\@split,\%aln);
		### set read equal to new
		$read_name = $split[0];		
	}
	
}

### you need to process last batch of alignments ###############################
collapse_alignments_same_ss(\%aln);
filter_and_print_alignments(\%aln);
print ">>> Finished Filtering Reads <<<\n";
################################################################################

### after you created the filtered SS_SAM file, we need to bam it ##############
print ">>> Converting Filtered SAM to BAM <<<\n";
!system("samtools view -bS $clean_sam > $bam") or die "[STDERR]: ";
system("rm $clean_sam");
print ">>> Finished Converting Filtered SAM to BAM <<<\n";
################################################################################

### after output file is converted to bam, sort it BY CHROM NAME ###############
print ">>> Sorting BAM <<<\n";
my $sort_prefix = $bam_name . "_sort";
system("samtools sort $bam $sort_prefix");
print ">>> Finished Sorting BAM <<<\n";
################################################################################

################################################################################	
### rename
my $bam_sorted = $bam_dir . $bam_name . "_sort" . $bam_ext;
system("mv $bam_sorted $bam");
##############################################################################

################################################################################
sub add_alignment {

	### reference to single alignment (in an array)
	my @alignment = @{shift(@_)};
	### reference to the hash of alignments
	my $alignments = shift;
	
	### increase the "index", which is stored in the hash. this index is simply used to index the alignments
	my $index = scalar(keys %{$alignments}) + 1;

	### split the alignment into mandatory (mand) and optional (opt) fields

	### put the mand fields in an array
	$alignments->{$index}{"mand"} = [@alignment[0..10]];

	### put the opt fields in a hash. eg. for NM:i:5, NM:i is the key and 5 is the value
	foreach (@alignment[11..$#alignment]) {
		my @split = split(':');
		$alignments->{$index}{"opt"}{$split[0].":".$split[1]} = $split[2];
	}
}

sub collapse_alignments_same_ss {
	
	### the alignments we need to sort through
	my $alignments = shift;

	### return unless there are at least 2 alignments
	return unless scalar(keys %{$alignments}) > 1;

	### we collapse alignments that have identical "flag","chromosome","position","cigar",and "opt" values (except SS:Z)
	my %UNIQUE;

	### loop through all alignments, create a "unique" key upon which we compress
	foreach my $index (sort {$a <=> $b} keys %{$alignments}) {
		
		### use the mandatory part of the read (cols 1 to 11) to determine uniqueness
		my $mand = join("\t",@{$alignments->{$index}{"mand"}});

		foreach my $opt_key (keys %{$alignments->{$index}{"opt"}}) {
			$UNIQUE{$mand}{$opt_key}{$alignments->{$index}{"opt"}{$opt_key}}++;			
		}
	}

	### return if the number of unique alignments is equivalent to the number of alignments
	return if scalar(keys %UNIQUE) == scalar(keys %{$alignments});

	### make a new hash containing the new unique alignments
	my %alignments_new;
	### add each unique read (indexed by "mand") into new hash
	foreach my $mand (keys %UNIQUE) {
		### @read_new will be the new unique read
		my @read_new = split('\t',$mand);
		
		foreach my $opt (keys %{$UNIQUE{$mand}}) {
			### for each optional argument, collapse all the values into one
			my @keys = sort keys %{$UNIQUE{$mand}{$opt}};
			### push the new values onto the end of the read 
			push(@read_new,$opt.':'.join(',',@keys));
		}
		add_alignment(\@read_new,\%alignments_new);
	}
	
	$alignments = \%alignments_new;
}

sub filter_and_print_alignments {	
	### This subroutine loops through all alignments and filters out reads that are greater than max strata and max_edit distance
	### We then print out alignments, appending 3 new tags (see below) and update the "mapqual" column to equal the # of total alignments this read has
	### OT: other alignments for this read
	### RS: relative rank of this strata (starting at 1)
	### NS: number of aln in this strata
	### TS: total number of strata
	my $alignments = shift;
	
	### first loop through all alignments and filter
	### if read has alignment, tally up the strata and create the "aln" string 
	my (%STRATA,%ALN_STRING);
	foreach my $index (sort {$a <=> $b} keys %{$alignments}) {
		### if read is unaligned, then just print it and we're done with this group of alignments.
		if ($alignments->{$index}{"mand"}->[1] =~ /^4$/) {
			print_alignment($alignments->{$index});
			return;
		} else {
			### if aln doesn't pass filter, delete alignment.
			if ($alignments->{$index}{"opt"}{"XA:i"} > $max_strata || $alignments->{$index}{"opt"}{"NM:i"} > $max_edit_distance) {
				delete $alignments->{$index};
			} else {
				### create the string that goes in the "OT" tag
				$ALN_STRING{$index} =	$alignments->{$index}{"mand"}->[2] . ':' .
							$alignments->{$index}{"mand"}->[3] . ':' .
							(($alignments->{$index}{"mand"}->[1] == 0 && '+') || '-') . ':' .
							$alignments->{$index}{"mand"}->[5] . ':' .
							$alignments->{$index}{"opt"}{"NM:i"};
				### add strata to tally
				$STRATA{$alignments->{$index}{"opt"}{"XA:i"}}++;
			}
		}
	}

	### sort the strata and calculate the relative order of strata (within the other strata)
	my (%STRATA_REL,$relative_order);
	foreach my $strata (sort {$a <=> $b} keys %STRATA) { $STRATA_REL{$strata} = ++$relative_order; }

	my $TOT_NUM_ALN = scalar(keys %{$alignments});
	my $TOT_STRATA = scalar(keys %STRATA_REL);
	foreach my $INDEX (sort {$a <=> $b} keys %{$alignments}) {
		my @ALN_STRING;
		foreach my $index (sort {$a <=> $b} keys %ALN_STRING) { push(@ALN_STRING,$ALN_STRING{$index}) if $index != $INDEX; }
		$alignments->{$INDEX}{"opt"}{"OT:Z"} = join(";",@ALN_STRING) if @ALN_STRING;
		$alignments->{$INDEX}{"mand"}->[4] = $TOT_NUM_ALN;
		$alignments->{$INDEX}{"opt"}{"NS:i"} = $STRATA{$alignments->{$INDEX}{"opt"}{"XA:i"}};
		$alignments->{$INDEX}{"opt"}{"RS:i"} = $STRATA_REL{$alignments->{$INDEX}{"opt"}{"XA:i"}};
		$alignments->{$INDEX}{"opt"}{"TS:i"} = $TOT_STRATA;
		print_alignment($alignments->{$INDEX});
	}
}

sub print_alignment {
	my $alignment = shift;

	print $CLEAN_SAM join("\t",@{$alignment->{"mand"}});
	foreach my $opt_field (sort keys %{$alignment->{"opt"}}) {
		print $CLEAN_SAM "\t", $opt_field, ":", $alignment->{"opt"}{$opt_field};
	}
	print $CLEAN_SAM "\n";
}
