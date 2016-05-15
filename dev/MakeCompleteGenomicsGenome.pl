#!/usr/bin/perl

use strict;
use Getopt::Long;

### help variable
my $help;
### variation file from complete genomics
my $varfile; ###"/home/jmtoung/Lab/rdd/complete_genomics/ASM_Build36/NA06985/GS06985-1100-36-ASM/GS00392-DNA_D01/ASM/var-GS06985-1100-36-ASM.tsv";
### output directory for genomes
my $outputdir; ### "/home/jmtoung/Lab/rdd/complete_genomics/ASM_Build36/NA06985/GS06985-1100-36-ASM/GS00392-DNA_D01/ASM";
### name of genomes
my $genome_name; ### "gm06985";
### reference genome
my $ref_genome; ### "/home/jmtoung/Lab/database/hg18/hg18.fa";
### list of chromosomes to output. if not specified, then all will be used.
my $chr_list; ### "chr19,chr17,chr22

my $result = GetOptions("help" => \$help, "varfile=s" => \$varfile, "outputdir=s" => \$outputdir, "genome_name=s" => \$genome_name, "ref_genome=s" => \$ref_genome, "chr_list=s" => \$chr_list);

### print help 
if (defined $help) {
	print "Usage: MakeCompleteGenomicsGenome.pl --varfile <var.tsv> --outputdir <outputdir> --genome_name <genome_name> --ref_genome <ref_genome.fa> --chr_list <chr1,chr2>\n\n";
	print "--varfile\tvariations file from Complete Genomics\n";
	print "--outputdir\tfolder where genome will be placed\n";
	print "--genome_name\tname of genomes and folder for genomes\n";
	print "--ref_genome\treference genome file\n";
	print "--chr_list\tchromosomes to be included (default is all)\n";
	exit;
}

### check if variables are undefined
unless(defined $varfile) {
	print "varfile is missing\n";
	exit;
} 
unless(defined $genome_name) {
	print "genome_name is missing\n";
	exit;
}
unless(defined $ref_genome) {
	print "ref_genome is missing\n";
	exit;
}
unless(-d $outputdir) {
	print "outputdir is missing (or doesn't exist)\n";
	exit;
}

### check if folders to contain genome indices already exist
chdir($outputdir);
if(-d $genome_name . "_1" || -d $genome_name . "_2") {
	print "folder for $genome_name exists already\n";
	exit;
} else {
	mkdir $genome_name . "_1";
	mkdir $genome_name . "_2";
}

### this chrom hash contains all info regarding chromosomes
my %chr;

### split $chr_list and put it in the %chrom hash
if (defined $chr_list) {
	my @chr_list = split(',',$chr_list);
	### for each chromosome specified in $chr_list, indicate that it's "use"
	foreach (@chr_list) {
		$chr{$_}{"use"}++;
	}
}

### get chromosomes from var file (column 4) 
{
	### awk command that gets the unique elements of column 4 from varfile
	my $awk = "awk '" . '$4~/chr[0-9XYM]+/ {print $4}' . "' $varfile | uniq |";
	open(AWK,$awk) or die "couldn't fork: $!\n";
	### index value keeps track of order of chromosomes in varfile
	my $index = 0;
	while(<AWK>) {
		chomp;
		$chr{$_}{"index"} = $index++;
		### if $chr_list was undef, it means all chromosomes are to be used
		unless (defined $chr_list) {
			$chr{$_}{"use"}++;
		}
	}
}

### check that chrY comes after chrX in the var file (due to the pseudo-auto region of chrY)
if (defined $chr{"chrY"} && $chr{"chrY"}{"index"} > $chr{"chrX"}{"index"}) {
	print "chrY > chrX\n";
	exit;
}

### open genome file, one for each allele
my %genome_fh;
open($genome_fh{"1"},">" . $genome_name . "_1/" . $genome_name . "_1.fa");
open($genome_fh{"2"},">" . $genome_name . "_2/" . $genome_name . "_2.fa");

### open map files, one for each allele
my %map_file;
open($map_file{"1"},">" . $genome_name . "_1/" . $genome_name . "_1.map");
open($map_file{"2"},">" . $genome_name . "_2/" . $genome_name . "_2.map");

### define sequence variables, one for each allele. this sequence variable keeps track of the sequence to be inserted
my %seq;
$seq{"1"} = "";
$seq{"2"} = "";

### define index variables, one for each allele. these index variables keep track of the length of the genome sequence
my %genome_length;

### this variable checks whether the header has ended
my $end_header;

### open and read varfile
open(VARFILE,$varfile);
while(<VARFILE>) {
	### if you get to the header row, set the $end_header variable to 1 to signal that we can process the rest of the file
	if (/>/) {
		$end_header = 1;
		next;
	}
	
	### if you're still in the header portion, skip
	next unless $end_header;

	### split the line
	chomp;
	my @line = split('\t');
	
	### get locus
	my $locus = $line[0];

	### get allele, if allele is all, change it to both 1 and 2
	my @allele;
	if ($line[2] eq "all") {
		@allele = ("1","2");
	} else {
		@allele = ($line[2]);
	}
	
	### get chromosome
	my $chrom = $line[3];
	unless (defined $chr{$chrom}{"use"}) {
		next;
	}
	
	### get start (add 1 to convert from 0-based coordinates to 1-based)
	my $start = $line[4] + 1;
	
	### get end
	my $end = $line[5];
	
	### get var type
	my $var_type = $line[6];
	
	### get ref
	my $ref = $line[7];
	
	### get alleleseq
	my $allele_seq = $line[8];
	
	### get total score
	my $score = $line[9] if defined $line[9];

	### get haplink
	my $haplink = $line[10] if defined $line[10];
	
	### get xref
	my $xref = $line[11] if defined $line[11];

	### check if this chrom has been seen, if not, print fasta header (">chr1")
	### also utilize this opportunity to flush out any remaining bases in the seq variable
	foreach my $allele (@allele) {
		unless (defined $chr{$chrom}{$allele}{"seen"}) {
			### print the leftovers (strings less than 50)
			print {$genome_fh{$allele}} $seq{$allele}, "\n" if (length($seq{$allele}) != 0);
			### clear seq
			$seq{$allele} = "";
			### clear genome_length
			$genome_length{$allele} = 1;
			### now print header
			print {$genome_fh{$allele}} ">$chrom", "\n";
			### indicate that this chromosome has been seen
			$chr{$chrom}{$allele}{"seen"}++;
		}
	}	

	### skip unless it's snp, del, ins, sub
	unless ($var_type eq "snp" || $var_type eq "del" || $var_type eq "ins" || $var_type eq "sub") {
		### pull sequence, do it here, not in foreach loop below b/c that would be repetitive
		my $ref_seq = pull_ref_seq($chrom,$start,$end);
		foreach my $allele (@allele) {
			### append sequence			
			$seq{$allele} .= $ref_seq;
			### print to map file
			print {$map_file{$allele}} $chrom, "\t", $start, "\t", $end, "\t", $genome_length{$allele}, "\t", "0", "\n";
			### update genome length
			$genome_length{$allele} += $end - $start + 1;
		}
	}
	
	if ($var_type eq "snp") {
		foreach my $allele (@allele) {
			### append snp allele to sequence (no need to pull out sequence from reference)
			$seq{$allele} .= $allele_seq;
			### print to map file
			print {$map_file{$allele}} $chrom, "\t", $start, "\t", $end, "\t", $genome_length{$allele}, "\t", "0", "\n";
			### update genome length (by 1)
			$genome_length{$allele}++;
		}
	} elsif ($var_type eq "del") {
		### since it's a deletion, no need to add to sequence. just print to map file
		foreach my $allele (@allele) {
			print {$map_file{$allele}} $chrom, "\t", $start, "\t", $end, "\t", "del", "\t", "0", "\n";
		}
	} elsif ($var_type eq "ins") {
		foreach my $allele (@allele) {
			### get insertion (cannot pull sequence from reference)
			$seq{$allele} .= $allele_seq;
			### get length of insertion
			my $length = length($allele_seq);
			### print to map file
			print {$map_file{$allele}} $chrom, "\t", $start, "\t", $end, "\t", "ins", "\t", $length, "\n";
			### update genome length to size of insertion
			$genome_length{$allele} += $length;
		}
	} elsif ($var_type eq "sub") {
		foreach my $allele (@allele) {
			### get sub (cannot pull sequence from reference)
			$seq{$allele} .= $allele_seq;
			### get length of sub
			my $length = length($allele_seq);
			### print to map file
			print {$map_file{$allele}} $chrom, "\t", $start, "\t", $end, "\t", "sub", "\t", $length, "\n";
			### update genome length to size of sub
			$genome_length{$allele} += $length;
		}
	
	}
	
	### print some of the sequence out (leaves 49 bases max in $seq{$allele})
	foreach my $allele (@allele) {
		$seq{$allele} = print_seq($seq{$allele},$genome_fh{$allele});
	}

}

### at the very end, print out remaining bases
my @allele = ("1","2");
foreach my $allele (@allele) {
	### print the leftovers (strings less than 50)
	print {$genome_fh{$allele}} $seq{$allele}, "\n" if (length($seq{$allele}) != 0);
}

### this function pulls the reference sequence out of the ref_genome using samtools
sub pull_ref_seq {
	my $chrom = shift(@_);
	my $start = shift(@_);
	my $end = shift(@_);
	
	### define region to be pulled
	my $region = $chrom . ":" . $start . "-" . $end;
	### pull region using samtools
	my $seq = `samtools faidx $ref_genome $region`;
	### split the output on new line character
	my @split = split('\n',$seq);
	### get rid of first line
	shift(@split);
	### return sequence concatenated
	return join("",@split);

}

### this function prints the seq provided, with 50 bases max on each line printed
sub print_seq {
	my $seq = shift(@_);
	my $fh = shift(@_);
	
	while(length $seq >= 50) {
		### pull first 50 bases of the string and substitute those bases with 'X'
		my $first_50 = substr($seq,0,50,'X');
		### print first 50 bases
		print {$fh} $first_50, "\n";
		### remove 'X' (essentially delete those 50 base pairs that were just printed
		$seq =~ s/X//g;
	}
	
	### return any left over string that was not printed (max length is 49!)
	return $seq;
}
