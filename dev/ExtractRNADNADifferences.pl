#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib "/ifs/apps/BioPerl-1.6.9/lib/perl5", '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use Database;
use File::Basename;
use PileupData;
use IlmnDNA;
use ComplementBase;
umask 0007;

my $home;
my $pileup;
my $region; ### added on 4/17/2012
my $tag;
my $dna;
my $strand_specific;
my $minqual;
my $alnDB;
### the following two options are soley so that we can say how many bases we analyzed for RDD
my $bed = 0;
my $min_total_reads = 10;

#############################################################################################

my $result = GetOptions(
	"home=s" => \$home,
	"pileup=s" => \$pileup,
	"region=s" => \$region,
	"tag=s" => \$tag,
	"dna=s" => \$dna,
	"strand_specific=i" => \$strand_specific,
	"minqual=i" => \$minqual,
	"alnDB=s" => \$alnDB,
	"bed=i" => \$bed,
	"min_total_reads=i" => \$min_total_reads
);

################################################################################
### This script outputs sites where RNA is different from DNA subject to minqual thresholds
### the 'bed' option means to just spit out the sites that were not skipped so we can tally up how many gigabases of sequenced was analyzed
################################################################################

(print STDERR "home:\t$home\n") && -d $home or die "[STDERR]: home not defined\n";
(print STDERR "pileup:\t$pileup\n") && -e $pileup or die "[STDERR]: pileup not defined\n";
(print STDERR "region:\t$region\n") if defined $region;
(print STDERR "tag:\t$tag\n") && defined $tag or die "[STDERR]: tag not defined\n";
(print STDERR "dna:\t$dna\n") && -e $dna or die "[STDERR]: dna not defined\n";
(print STDERR "strand_specific:\t$strand_specific\n") && defined $strand_specific or die "[STDERR]: strand specific not defined\n";
(print STDERR "minqual:\t$minqual\n") && defined $minqual or die "[STDERR]: minqual not defined\n";
(print STDERR "alnDB:\t$alnDB\n") && -d $alnDB or die "[STDERR]: alnDB not defined\n";
(print STDERR "bed:\t$bed\n") or die "[STDERR]: bed not defined\n";
(print STDERR "min_total_reads:\t$min_total_reads\n") or die "[STDERR]: min_total_reads not defined\n";

my ($regionPart,$regionStart,$regionEnd);
if (defined $region) {
	($regionPart,$regionStart,$regionEnd) = split('-',$region);
}

my ($pileup_name,$pileup_dir,$pileup_ext) = fileparse($pileup,'\.pileup(\.gz)?');
my $PILEUP;
if (defined $regionStart) {
	$PILEUP = PileupData->new($pileup,$minqual,$regionStart);
} else {
	$PILEUP = PileupData->new($pileup,$minqual);
}

### LOOK UP bamID ##############################################################
my $bam = $pileup_dir;
$bam =~ s/(\.unique)?\.pileup(-noBAQ)?\//\.bam/;
my $home_regex = $home;
$home_regex =~ s/\//\\\//g;
$bam =~ s/$home_regex//g;
$bam =~ s/^\///;
my $alnID = Database->new($alnDB)->lookup(1,$bam,0);
defined $alnID or die "[STDERR]: alnID for '$bam' not defined in '$alnDB'\n";
################################################################################

### LOAD DNA ###################################################################
my $DNA = IlmnDNA->new($dna);
################################################################################

my @strands;
if ($strand_specific) { push(@strands,['+'],['-']); }
else { push(@strands,['+','-']); }

### MAKE OUTPUT ################################################################
unless ($bed) {
	my $output_file;
	if (defined $region) {
		my $output_dir = $pileup_dir = $pileup_dir . $pileup_name . ".$tag";
		system("mkdir $output_dir") unless -d $output_dir;
		$output_file = $output_dir . "/" . $pileup_name . ".$tag" . ".$regionPart" . ".txt";
	} else {
		$output_file = $pileup_dir . $pileup_name . ".$tag.txt";
	}	
	print STDERR $output_file, "\n";
	open(OUTPUT,">$output_file") or die "[STDERR]: can't open $output_file: $!\n";
	my $old_fh = select(OUTPUT); $|++; select($old_fh);
}

my %SKIP; ### hash of info skipped sites

### AN ARRAY OF THE DNA SITES THAT OVERLAP CURRENT SITE #######################
my @DNA_overlap;
################################################################################
while($PILEUP->has_next) {
	my $chrom = $PILEUP->get_chrom; 
	my $position = $PILEUP->get_position;

	if (defined $region) {
		last if $position > $regionEnd;
	}
	
	my $ref_base = $PILEUP->get_ref_base;

	goto NEXT if $ref_base eq 'N';

	### if the start is equal to or less than the current site, push it to @DNA_overlap
	while($DNA->has_next && $DNA->get_start <= $position) { 
###		print STDERR "."; ### turned this off on 5/2/2012
		push(@DNA_overlap,$DNA->get_next_line) if $DNA->get_next_line->[2] >= $position; ### added the if clause on 4/17
		$DNA->update_next_line;
	}

	### for each thing in @DNA_overlap check that the end is equal to or greater than current site
	my @DNA_overlap_temp;
	foreach my $DNA_overlap (@DNA_overlap) {
		push(@DNA_overlap_temp,$DNA_overlap) if $DNA_overlap->[2] >= $position;
	}
	@DNA_overlap = (); @DNA_overlap = @DNA_overlap_temp;

	### remove 'dna' if there's more than one polymorphism
	if (@DNA_overlap > 1) {
		my @DNA_overlap_rm_dna;
		foreach my $DNA_overlap (@DNA_overlap) {
			push(@DNA_overlap_rm_dna,$DNA_overlap) unless $DNA_overlap->[3] eq 'dna';
		}
		@DNA_overlap = @DNA_overlap_rm_dna;
	}

	### only want HOMOZYGOUS (concordant) genotypes, skip everything else!
	### WANT: snp: homozygous_concordant, indel: homozygous/reference, deletion: homozygous/reference, dna
	my ($var_type,$genotype,$DNA_base);
	my $SKIP = 0; ### to determine if we want to skip
	if (@DNA_overlap > 1) {
		$SKIP = 1; ### skip if there's more than one polymorphism here
		my (@var_type,@genotype,@DNA_base);
		foreach my $DNA_overlap (sort {$a->[3] cmp $b->[3]} @DNA_overlap) {
			my $var_type = $DNA_overlap->[3];
			push(@var_type,$var_type);
			push(@genotype,get_genotype($var_type,$DNA_overlap));
			push(@DNA_base,get_DNA_base($var_type,$DNA_overlap,$ref_base));	
		}
		$var_type = join(";",@var_type);
		$genotype = join(";",@genotype);
		$DNA_base = join(";",@DNA_base);
	} elsif (@DNA_overlap == 1 && $DNA_overlap[0]->[3] eq 'dna') {
		$var_type = 'dna';
		$genotype = 'dna';
		$DNA_base = $ref_base;
	} elsif (@DNA_overlap == 1) {
		$var_type = $DNA_overlap[0]->[3];
		$genotype = get_genotype($var_type,$DNA_overlap[0]);
		$DNA_base = get_DNA_base($var_type,$DNA_overlap[0],$ref_base);
		$SKIP = 1 if $var_type =~ /snp/ && $genotype ne 'homozygous_concordant';
		$SKIP = 1 if $genotype =~ /heterozygous/ && ($var_type =~ /insertion/ || $var_type =~ /deletion/); 
	} elsif (@DNA_overlap == 0) {
		$var_type = 'monomorphic';
		$genotype = 'monomorphic';
		$DNA_base = $ref_base;
	}

	my $has_rdd = 0;	

	goto SKIP if $SKIP;
	
	### get the counts of the bases
	my $DNA_base_counts;
	if (@DNA_overlap == 1) {
		my $length = @{$DNA_overlap[0]} - 1;
		$DNA_base_counts = join("\t", @{$DNA_overlap[0]}[9..$length]);
	}

	foreach my $strands (@strands) {
		my $bases = $PILEUP->get_bases($strands);
		
		my @total_bases;
		foreach my $base (@{$bases}) {
			push(@total_bases,$base) unless $base eq 'S';
		}
		my $total_count = $PILEUP->get_base_count($strands,\@total_bases);
		
		if ($bed) {
			print $chrom, "\t", $position - 1, "\t", $position, "\t", $alnID, "\n" if $total_count >= $min_total_reads;
			goto NEXT;
		}
		
		foreach my $base (@{$bases}) {
			next if $base eq 'S';
			
			my ($REF_BASE,$DNA_BASE);
			if ($strands->[0] eq '-') { 
				$REF_BASE = complement_base($ref_base);
				$DNA_BASE = complement_base($DNA_base);
			} else { 
				$REF_BASE = $ref_base; 
				$DNA_BASE = $DNA_base;
			}

			next if $base eq $DNA_BASE;
			
			my $base_count = $PILEUP->get_base_count($strands,[$base]);

			my @OUTPUT;
			push(@OUTPUT,$chrom,$position,join(",",@{$strands}),$DNA_BASE,$total_count,$base,$base_count,$var_type,$genotype,$REF_BASE,$alnID);
			push(@OUTPUT,$DNA_base_counts) if defined $DNA_base_counts;
			print OUTPUT join("\t",@OUTPUT), "\n";
			$has_rdd++;
		}
	}

	SKIP:
	$SKIP{$SKIP}{$has_rdd}{$var_type}{$genotype}{$DNA_base}++;
	
	NEXT: 

	$PILEUP->parse_next_line;
}
close(OUTPUT);

### exit if $bed; ### commented this out on 4/17

foreach my $SKIP (sort {$a<=>$b} keys %SKIP) {
	foreach my $has_rdd (sort {$a<=>$b} keys %{$SKIP{$SKIP}}) {
		foreach my $var_type (sort keys %{$SKIP{$SKIP}{$has_rdd}}) {
			foreach my $genotype (sort keys %{$SKIP{$SKIP}{$has_rdd}{$var_type}}) {
				foreach my $DNA_base (sort keys %{$SKIP{$SKIP}{$has_rdd}{$var_type}{$genotype}}) {
					print STDERR "tally\t", $pileup, "\t", $tag, "\t", $alnID, "\t", $SKIP, "\t", $has_rdd, "\t", $var_type, "\t", $genotype, "\t", $DNA_base, "\t", $SKIP{$SKIP}{$has_rdd}{$var_type}{$genotype}{$DNA_base}, "\n";
				}
			}
		}
	}
}

sub get_genotype {
	my ($var_type,$DNA_overlap) = @_;
	
	if ($var_type =~ /deletion/ || $var_type =~ /insertion/) {
		return $DNA_overlap->[7]; ###  (homozygous,heterozygous,reference)
	} elsif ($var_type =~ /snp/) {
		return $DNA_overlap->[7] . "_" . $DNA_overlap->[8]; ### combination of calls from two algorithm parameters
	}
}

sub get_DNA_base {
	my ($var_type,$DNA_overlap,$ref_base) = @_;
	
	if ($var_type =~ /deletion/) {
		my $genotype = get_genotype($var_type,$DNA_overlap);
		if ($genotype eq 'reference') {
			return $ref_base;
		} elsif ($genotype eq 'homozygous') {
			return 'X';
		} elsif ($genotype eq 'heterozygous') {
			return $ref_base . '|' . 'X';
		}
	} elsif ($var_type =~ /insertion/) {
		my $genotype = get_genotype($var_type,$DNA_overlap);
		if ($genotype eq 'reference') {
			return $ref_base;
		} elsif ($genotype eq 'homozygous') {
			return $ref_base . $DNA_overlap->[5];
		} elsif ($genotype eq 'heterozygous') {
			return $ref_base . '|' . $ref_base . $DNA_overlap->[5];
		}
	} elsif ($var_type =~ /snp/) {
		my $genotype = get_genotype($var_type,$DNA_overlap);
		if ($genotype eq 'homozygous_concordant') {
			my @alleles = split('',$DNA_overlap->[5]);
			$alleles[0] eq $alleles[1] or die "[STDERR]: something wrong with homozygous concordant alleles...\n";
			return $alleles[0];
		}
		return $DNA_overlap->[5] if $DNA_overlap->[5] eq $DNA_overlap->[6]; ### if calls are concordant just return one
		return $DNA_overlap->[5] . "," . $DNA_overlap->[6]; ### if calls are discordant
	}
}
