#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

### bam file
my $bam;
my $ss_map;

my $result = GetOptions("bam=s" => \$bam, "ss_map=s" => \$ss_map);

### check that $bam ends in '.bam' and is an absolute file name ################
my ($bam_name, $bam_dir, $bam_ext) = fileparse($bam,'\.bam');
$bam_ext eq '.bam' or die "[STDERR]: '$bam' does not end in '.bam'\n";
substr($bam_dir,0,1) eq '/' or die "[STDERR]: '$bam' is not absolute file name\n";
chdir($bam_dir) or die "[STDERR]: cannot chdir '$bam_dir': $!\n";
### open the bam file to read using `samtools view`
my $BAM = "samtools view -h " . $bam . " |";
open(BAM,$BAM) or die "[STDERR]: couldn't fork $BAM: $!\n"; 
### open the file to write sam to
my $SAM = $bam_dir . $bam_name . "_conv.sam";
open(SAM,">".$SAM) or die "[STDERR]: cannot open '$SAM': $!\n";
################################################################################

### store map file in hash #####################################################
my %ss_map;
open(SS_MAP,$ss_map) or die "[STDERR]: cannot open '$ss_map': $!\n";
while(<SS_MAP>) {	
	### if it's the header, skip
	next if (/^###/);
	chomp;
	my @split = split('\t');
	### $split[0] is the ss_id and $split[1] is the splice structure
	$ss_map{$split[0]} = $split[1];
}
################################################################################

print ">>> Converting Coordinates <<<\n";
### convert coordinates in the bam file to hg18 coordinates (from the splice library coordinates)
while(<BAM>) {
	chomp;

	die "[STDERR]: $_\n" if $_ eq "[bam_header_read] EOF marker is absent. The input is probably truncated.";
	
	my @split = split('\t');

	### if it's a header && the @SQ is for a splice read (we don't want them anymore!!!), then skip
	if ($split[0] eq '@SQ' && $split[1] =~ /SN:ss[0-9]+/) {
		next;
	### if it's the header line OR if this read has no alignments (flag is 4) OR if the seq_id is not an integer just print and skip
	} elsif (/^@/ || $split[1] =~ /^4$/ || !($split[2] =~ /^ss[0-9]+$/)) {
		print SAM $_, "\n";
	} else {
	### otherwise, we have to change the seq_id, the position, the cigar string, and append the old seq_id to the end
		### append the ss_id to the end as an optional field
		push(@split,"SS:Z:".$split[2]);
		### change ss_id ($split[2]) to chrom, position ($split[3]) to hg18 coordinates, and cigar ($split[5])
		($split[2],$split[3],$split[5]) = convert_position($split[2],$split[3],\%ss_map,length($split[9]));
		### print the new alignment (@split), but also check for undef
		if (!defined $split[2] || !defined $split[3] || !defined $split[5]) {
			print STDERR "[STDERR]: undef_values\t$_\n";
		} else {
			print SAM join("\t",@split), "\n";
		}
	}
}
print ">>> Finished Converting Coordinates <<<\n";
close(SAM);
close(BAM);

### after you created the SS_SAM file, we need to convert that to bam.
print ">>> Converting SAM to BAM <<<\n";
system("samtools view -bS -o $bam $SAM");
print ">>> Finished Converting SAM to BAM <<<\n";

### delete the sam file
print ">>> Deleting SAM <<<\n";
system("rm $SAM");
print ">>> Finished Deleting SAM <<<\n";

### after output file is converted to bam, sort it BY READ NAME
print ">>> Sorting BAM by Read Name <<<\n";
my $sort_prefix = $bam_name . "_sortn";
system("samtools sort -n $bam $sort_prefix");
print ">>> Finished Sorting BAM by Read Name <<<\n";

### rename
my $BAM_sorted = $bam_dir . $bam_name . "_sortn.bam"; 
system("mv $BAM_sorted $bam");

################################################################################

sub convert_position {
	
	### this is the ss_id of the splice_site (e.g. ss1, ss2, ss3, etc.)
	my $ss_id = shift(@_);
	### this is the starting position of the alignment to this splice_read (as determined from bowtie/bwa, whatever)
	my $align_position = shift(@_);
	### this is the hash that contains the map between seq_id and splice structure
	my $ss_map = shift(@_);
	### length of this read
	my $read_length = shift(@_);
	
	### get the structure of the splice junction from the ss_map hash
	my $ss_structure = $$ss_map{$ss_id};

	########################################################################
	### split the structure to get (1) $chrom and (2) $positions, which will then be further split to get @starts and @ends
	### e.g. of a structure: chr19|45232038,45232086;45464754,46465895|45232254,45232302
	my ($chrom, @positions) = split('\|',$ss_structure);
	my (@starts, @ends);
	### for both the left and right side, split parts and add to @starts and @ends array
	foreach my $side (@positions) {
		### split the structure into parts (sometimes there will be multiple parts of the left/right side
		my @parts = split(';',$side);
		### for each start_end pair of a portion of the splice structure, push it to a starts and ends array
		foreach my $part (@parts) {
			### split to get start and end
			my ($start,$end) = split(',',$part);
			push(@starts,$start);
			push(@ends,$end);
		}
		
	}
	########################################################################

	### this is the converted align position we want to return
	my $converted_align_position;
	### this is the cigar string (that we want to return!)
	my $cigar;
	### there are two loops; the first loop ($i) is for finding where the starting position begins
	### the second loop is for creating the cigar string and it begins once we find where the starting position begins!
	for (my $i = 0; $i <= $#starts; $i++) {
		### this is the length of the entire exon we are currently interrogating
		my $exon_length = $ends[$i] - $starts[$i] + 1;
		
		### if the align_position is less than the length of the entire exon, then this means that the starting position of the alignment must be within this exon.
		### if not, we deduct the length of this entire exon from the $align_position and move onto the next exon.
		if ($align_position <= $exon_length) {
			### calculate the $converted_align_position.
			### $start_piece is a variable that we'll use to calculate the length of this "piece" we're working on. normally, $start_piece should just be equal to 
			### the start of the exon ($starts[$i]), but for the exon to which the left most part of the read aligns, it is not
			my $start_piece = $converted_align_position = $starts[$i] + $align_position - 1;

			### LOOP 2	
			### $read_length starts out as the length of the read, but will gradually decrease as we match the read up to pieces of the exon.
			### once $read_length is 0, we're done.
			while ($read_length > 0) {
				### calculate the maximal length of this exon. for the first exon to which the read aligns, the start is the start of the alignment, 					
				### not the start of the exon...(explained above).
				my $aligned_length = $ends[$i] - $start_piece + 1;
				### we actually need to compare it to the $read_length (or the remaining part of the read that isn't "aligned").
				### if $read_length is greater, we need to go onto next exon. if not we're done
				if ($read_length > $aligned_length) {
					### for the cigar string, add the number that matched this exon in addition to the gap between this exon and the next
					$cigar .= $aligned_length . "M" . ($starts[$i+1]-$ends[$i]-1) . "N";
					### deduct from the $read_length the number that has been aligned
					$read_length -= $aligned_length;
				} else {
					$cigar .= $read_length . "M";
					return ($chrom, $converted_align_position, collapse_cigar($cigar));	
				}
				### if we have to go to next exon, be sure to increase $i and update the $start_piece.
				$start_piece = $starts[$i+1];
				$i++;
			}
		} else {
			$align_position -= $exon_length;
		}
	}
}

sub collapse_cigar {
	### collapses cigar string by joining adjacent M's and N's
	
	my $cigar = shift(@_);
	### split the cigar string by M/N but keep the delimiter
	my @split = split(/(M|N)/,$cigar);

	### keeps track of the nums. shift to take first one
	my $num = shift(@split);
	### keeps track of last type (M or N). shift to take first one.
	my $type = shift(@split);
	### new cigar string to return	
	my @new_cigar;
	
	### pull two pieces from the @split each time (a number and a type)
	while (@split) {
		### compare compare_num to num
		my $compare_num = shift(@split);
		### compare compare_type to type
		my $compare_type = shift(@split);
		
		### if compare_type equals type (the type of the last pair), then...
		if ($compare_type eq $type) {
			### add to num and continue
			$num += $compare_num;
		} else {
			### if not, push the accumulated num into new cigar
			push(@new_cigar,$num);
			### push type into new cigar
			push(@new_cigar,$type);
			### set type to new type
			$type = $compare_type;
			### set num to new num
			$num = $compare_num;
		}
	}
	
	### push $num and $type to @new cigar b/c exited before could do so for last pair
	push(@new_cigar,$num);
	push(@new_cigar,$type);

	return join("",@new_cigar);
}
