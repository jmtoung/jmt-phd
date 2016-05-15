#!/usr/bin/perl -w

use strict;
use POSIX;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev', '/gpfs/fs121/h/toung/oldhome/dev';
use File::Basename;
use PullGenomeSequence;
use Bio::SeqIO;
use Bio::Seq;
use ComplementBase;

my $HOME;
my $blat;
my $files;
my $index;
my $tag;
my $flanking_distance;
my $kmer = 12;
my $minIdentity;
my $repMatch = 2253;
my $maxBlocks = 3;
my $gapPen = 1;

my $options = GetOptions(
	"home=s" => \$HOME,
	"blat=s" => \$blat,
	"files=s" => \$files,
	"index=s" => \$index,
	"tag=s" => \$tag,
	"flanking_distance=s" => \$flanking_distance,
	"kmer=i" => \$kmer,
	"minIdentity=i" => \$minIdentity,
	"repMatch=i" => \$repMatch,
	"maxBlocks=i" => \$maxBlocks,
	"gapPen=i" => \$gapPen
);

(print STDERR "HOME:\t$HOME\n") && -d $HOME or die "[STDERR]: undefined HOME\n";
(print STDERR "blat:\t$blat\n") && -e $blat or die "[STDERR]: cant find blat\n";

my @files = split(',',$files);
$files = join(" ",@files);

foreach my $file (@files) { 
	(print STDERR "file:\t$file\n") && -e $file or die "[STDERR]: $file doesn't exist: $!\n";
}

(print STDERR "index:\t$index\n") && defined $index or die "[STDERR]: undefined $index\n";
(print STDERR "tag:\t$tag\n") && defined $tag or die "[STDERR]: undefined $tag\n";
(print STDERR "kmer:\t$kmer\n") && defined $kmer or die "[STDERR]: undefined $kmer\n";
(print STDERR "minIdentity:\t$minIdentity\n") && defined $minIdentity or die "[STDERR]: undefined $minIdentity\n";
(print STDERR "repMatch:\t$repMatch\n") && defined $repMatch or die "[STDERR]: undefined $repMatch\n";
(print STDERR "maxBlocks:\t$maxBlocks\n") && defined $maxBlocks or die "[STDERR]: undefined $maxBlocks\n";
(print STDERR "gapPen:\t$gapPen\n") && defined $gapPen or die "[STDERR]: undefined $gapPen\n";

my @flanking_distance = split(',',$flanking_distance);

foreach my $fl_dist (@flanking_distance) { 
	(print STDERR "flanking_distance:\t$fl_dist\n") && $fl_dist =~ /[0-9]+/ or die "[STDERR]: flanking distance invalid\n"; 
}

### make fasta files with sequences for each flanking distance
my %fasta;
foreach my $fl_dist (@flanking_distance) {
	my $fasta_file = "BlatRDDSeqExplainRDD-$tag-$fl_dist.fasta";
	print STDERR "fasta_file:\t$fasta_file\n";
	$fasta{$fl_dist}{'fasta_file'} = $fasta_file;
	$fasta{$fl_dist}{'fasta_obj'} = Bio::SeqIO->new(-file => ">$fasta{$fl_dist}{'fasta_file'}", -format => 'fasta');
}

### read the unique sites from the files and create sequences to blat
my $command = "cat $files | sort +0 -1 +1n -2 +2 -3 +3 -4 +5 -6 | uniq |";
print STDERR "cat_command:\t$command\n";

open(COMMAND,$command) or die "[STDERR]: can't fork $command: $!\n";
my %current; ### hold current RDD
while(<COMMAND>) {
	chomp;
	my @split = split('\t');
	$split[2] eq '+,-' or die "[STDERR]: strand doesn't equal '+,-'\n";
	
	if (is_equal(\@split,\%current)) { 
		add_to_current(\@split,\%current);
	} else {
		make_sequences(\%current,\@flanking_distance,\%fasta,$index);
		%current = ();
		add_to_current(\@split,\%current);
	}
}
make_sequences(\%current,\@flanking_distance,\%fasta,$index);

## RUN BLAT ON EACH OF THE FASTA FILES ########################################
my %blat_results;
foreach my $fl_dist (@flanking_distance) {
	my $blat_output = $fasta{$fl_dist}{'blat'} = "BlatRDDSeqExplainRDD-$tag-$fl_dist.pslx";

	my $blat_command = "time $blat -stepSize=5 -minIdentity=$minIdentity -repMatch=$repMatch -noHead -out=pslx $index $fasta{$fl_dist}{'fasta_file'} $blat_output";
	print STDERR $blat_command, "\n";
	system($blat_command);

	### parse results from blat output
	open(BLAT,$blat_output) or die "[STDERR]: can't open $blat_output: $!\n";

	while(<BLAT>) {
		chomp;
		my ($matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts,$qSeq,$tSeq) = split('\t');
		$qSeq =~ s/,$//; $tSeq =~ s/,$//; $blockSizes =~ s/,$//; $qStarts =~ s/,$//; $tStarts =~ s/,$//;

		### query info
		my ($chrom,$position,$flank_dist,$ref_base,$rdd_base,$mol) = split(';',$qName);
		my $region_start = $position - $flank_dist - 1; ### make it 0-based so can compare to $tStart!  ### the start of the region we blatted
			
		### next if the start of the target region is the same as the start of the region we blatted
		next if ($chrom eq $tName && (($tStart - $qStart) == $region_start));

		### calculate the number of mismatches (M) allowed. M = (readlength + 2)/kmer -2 + 2 (+2 to account for clipping of GSNAP)
		my $misMatches_allowed = floor(($qSize + 2)/$kmer);

		### step 1, don't allow for insertions in the query
		next if $qNumInsert > 0;
		
		### step 2, require block count be no greater than $maxBlocks
		next if $blockCount > $maxBlocks;

		## 0-based RDD position on the read, forward strand, the same as $flankDist
		my $qRDDPos = $flank_dist;
		## qBlockStart and qBlockEnd are with regard to the reverse strand, make qRDDPos the same way
		$qRDDPos = $qSize - 1 - $qRDDPos if $strand eq '-';
		
		## figure out which block contains the RDD site
		my @blockSizes = split(",", $blockSizes);
		my @qStarts = split(",", $qStarts);
		my @tStarts = split(",", $tStarts);
		my @tSeq = split(",", $tSeq);

		my $rddBlock = -1;
		for (my $i = 0; $i < $blockCount; $i++) {
			## blockStarts coordinates are with regard to reverse strand if query is mapped to reverse strand
			my $qStart = $qStarts[$i];
			my $qEnd = $qStart + $blockSizes[$i] - 1;

			## if RDD position is located within the current block, check to see if the aligned base in target is the same as RNA allele
			if ($qRDDPos >= $qStart && $qRDDPos <= $qEnd) {
				$rddBlock = $i;
				my $blockRDDPos = $qRDDPos - $qStarts[$i];
				my $tRDDBase = substr($tSeq[$i], $blockRDDPos, 1);
				$tRDDBase = complement_base($tRDDBase) if ($strand eq "-");
				$tRDDBase = uc($tRDDBase);
		    	
		    		my $rdd_explained = 0;
				$rdd_explained = 1 if $tRDDBase eq $rdd_base;
				
				
				my $num_misMatches = $qSize - ($matches + $repMatches) + $gapPen * $tNumInsert;
	
				if ($num_misMatches <= $misMatches_allowed) {

					$blat_results{$chrom}{$position}{$ref_base}{$rdd_base}{$fl_dist}{$rdd_explained}{$num_misMatches}++;
				}
			}
		}
	}	
	close(BLAT);
}

### open each file
foreach my $file (@files) {
	open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";
	my ($file_name,$file_dir,$file_ext) = fileparse($file,'\.txt');
	my $output = $file_dir . $file_name . "_blat-explRDD" . $file_ext;
	open(OUTPUT,">".$output) or die "[STDERR]: can't open $output: $!\n";
	
	while(<FILE>) {
		chomp;
		my @split = split('\t');
		### had to change this for ~/Lab/work/2011-12-25_MergeRNADNADifferences-noBAQ/twins_12004_12750/keep-output/min_rdd_level_5pct_single_nuc_changes/blat-hg18
		### for example, position is actually 2, etc.
		my $chrom = $split[0];
		my $position = $split[1];
		my $ref_base = $split[3];
		my $rdd_base = $split[5];

		foreach my $fl_dist (@flanking_distance) {
			foreach my $rdd_explained (0,1) {
				my @num_misMatches;
				foreach my $num_misMatch (sort {$a<=>$b} keys %{$blat_results{$chrom}{$position}{$ref_base}{$rdd_base}{$fl_dist}{$rdd_explained}}) {
					push(@num_misMatches,$num_misMatch . ':' . $blat_results{$chrom}{$position}{$ref_base}{$rdd_base}{$fl_dist}{$rdd_explained}{$num_misMatch});
				}
				push(@num_misMatches,'-') unless @num_misMatches;
				push(@split,join(";",@num_misMatches));
			}
		}
		print OUTPUT join("\t",@split), "\n";
	}
	close(FILE);
	close(OUTPUT);
}

sub is_equal {
	my ($split,$current) = @_;
	
	### if hash is empty (like for first one), return equal
	return 1 if !exists $current->{'chrom'};
	
	return 0 unless $current->{'chrom'} eq $split->[0];
	return 0 unless $current->{'position'} eq $split->[1];
	return 0 unless $current->{'strands'} eq $split->[2];
	return 0 unless $current->{'ref_base'} eq $split->[3];
	
	return 1;
}

sub add_to_current {
	my ($split,$current) = @_;
	
	$current->{'chrom'} = $split->[0];
	$current->{'position'} = $split->[1];
	$current->{'strands'} = $split->[2];
	$current->{'ref_base'} = $split->[3];
	$current->{'rdd_base'}{$split->[5]}++;
}

sub make_sequences {
	my ($current,$flank_dist,$fasta,$index) = @_;
	
	foreach my $fl_dist (@{$flank_dist}) {
		my $display_id = join(";",$current->{'chrom'},$current->{'position'},$fl_dist,$current->{'ref_base'});
		my $region = $current->{'chrom'} . ":" . ($current->{'position'} - $fl_dist) . "-" . ($current->{'position'} + $fl_dist);
		my $sequence = pull_sequence($region,$index);
		### sequence may not be defined if at ends of chromosome
		next unless $sequence && length($sequence) == ($fl_dist*2 + 1);

		my @rdd_bases = sort keys %{$current->{'rdd_base'}};
		my $rdd_bases = join(",",@rdd_bases);
		
		my @sequence = split('',$sequence);
		foreach my $rdd_base (@rdd_bases) {
			my @RNA_sequence = @sequence;
			if ($rdd_base eq 'X') {
				$RNA_sequence[$fl_dist] = '';	
			} else {
				$RNA_sequence[$fl_dist] = $rdd_base;
			}
			my $RNA_sequence = join('',@RNA_sequence);
			my $rna_obj = Bio::Seq->new(-seq => $RNA_sequence, -display_id => $display_id . ";$rdd_base;RNA", -desc => "", -alphabet => "dna");
			$fasta->{$fl_dist}{'fasta_obj'}->write_seq($rna_obj);
		}
	}
}
