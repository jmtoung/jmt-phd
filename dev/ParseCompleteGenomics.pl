#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use Cwd;

### var file
my $cg_file; # = "/home/jmtoung/Lab/rdd/complete_genomics/Diversity/ASM_Build36_1.10.0/NA12004/GS12004-1100-36-ASM/GS00392-DNA_B02/ASM/var-GS12004-1100-36-ASM.tsv";

my $result = GetOptions("cg_file=s" => \$cg_file);

################################################################################
### This script parses a var file from Complete Genomics into output ready for RDD calling
################################################################################

################################################################################
my ($cg_name, $cg_dir, $cg_ext) = fileparse($cg_file,'\.[a-zA-Z]+$');
################################################################################

################################################################################
open(CG_FILE,$cg_file) or die "$cg_file doesn't exist: $!\n";
if ($cg_dir eq './') { $cg_dir = getcwd . '/'; }
my %HASH; ### hash to store data
my %DIR;
foreach my $type ('snp','ins','del','sub') {
	$DIR{$type} = $cg_name . '_' . $type;
	chdir($cg_dir);
	unless(-d $DIR{$type}) { mkdir $DIR{$type}; }
}
################################################################################

################################################################################
### loop through file and store any data related to 'snp','ins','del','sub' in %HASH
my $header_passed = 0;
while(<CG_FILE>) {


	### skip through header
	unless ($header_passed) {
		if (/^>/) { $header_passed = 1; next; }
		else { next; }
	}
	
	chomp;	
	my @line = split('\t');
	next if $line[2] eq 'all'; ### seems that all variations described above ('snp','ins','del','sub') do not come with 'all' as the allele.
	
	if ($line[6] eq 'snp') {
		### $HASH->'snp'->'chrom1'->'position'->'allele' = 'genotype'
		$HASH{'snp'}{$line[3]}{$line[5]}{$line[2]} = $line[8];
	} elsif ($line[6] eq 'ins') {
		$HASH{'ins'}{$line[3]}{$line[5]}{$line[2]} = $line[8];
	} elsif ($line[6] eq 'del') {
		### for every base position that's deleted, give it a 'D' allele meaning deletion
		### note the data is stored in a diff way from 'snp' and 'ins'
		foreach (my $position = ($line[4] + 1); $position <= $line[5]; $position++) { 
			$HASH{'del'}{$line[3]}{$position}++; 
		}
	} elsif ($line[6] eq 'sub') {
		### if the length of the substitution is the same as that reference, then we can map each base 1:1
		if (length($line[7]) == length($line[8])) {
			my @ref = split('',$line[7]); ### the reference sequence that is replaced
			my @var = split('',$line[8]); ### the replacement sequence
			foreach (my $position = ($line[4] + 1); $position <= $line[5]; $position++) {
				if ($ref[0] ne $var[0]) { ### only update if the ref and var are different 
					$HASH{'sub'}{$line[3]}{$position}{$line[2]} = shift(@var);
				}
			}
		### if the length of the ref and sub are different, then for every position, just put the entire sub sequence there
		} else {
			foreach (my $position = ($line[4] + 1); $position <= $line[5]; $position++) {
				$HASH{'sub'}{$line[3]}{$position}{$line[2]} = $line[8];
			}
		}
	}
}
################################################################################

################################################################################
### now loop through each type of variation 'snp','ins', etc and print them out.
foreach my $type (keys %HASH) {
	chdir($cg_dir);
	chdir($DIR{$type});
	foreach my $chrom (keys %{$HASH{$type}}) {
		open(OUTPUT,'>'.$cg_name.'_'.$type.'_'.$chrom.'.txt');
		select(OUTPUT);
		foreach my $position (sort {$a <=> $b} keys %{$HASH{$type}{$chrom}}) {
			### chrom	position	type
			print $chrom, "\t", $position, "\t", $type, "\t";
			### if it's a snp, need to print both alleles together
			### however, both alleles are not necessarily defined, so need to check for existance
			if ($type eq 'snp') {
				if (defined $HASH{$type}{$chrom}{$position}{1}) { print $HASH{$type}{$chrom}{$position}{1}; } else { print "*"; }
				if (defined $HASH{$type}{$chrom}{$position}{2}) { print $HASH{$type}{$chrom}{$position}{2}, "\n"; } else { print "*\n"; }
			### if it's an insertion, print all (unique) insertions across the 2 alleles separated by commas. 
			} elsif ($type eq 'ins') {
				my %insertion;
				if (defined $HASH{$type}{$chrom}{$position}{1}) { $insertion{$HASH{$type}{$chrom}{$position}{1}}++; }
				if (defined $HASH{$type}{$chrom}{$position}{2}) { $insertion{$HASH{$type}{$chrom}{$position}{2}}++; }
				print join(',',keys %insertion), "\n";
			### if it's a deletion, just print 'D' for deletion
			} elsif ($type eq 'del') {
				print "D\n";
			### if it's a substitution, just like insertion, print all (unique) substitutions across the 2 alleles separated by commas
			} elsif ($type eq 'sub') {
				my %substitution;
				if (defined $HASH{$type}{$chrom}{$position}{1}) { $substitution{$HASH{$type}{$chrom}{$position}{1}}++; } 
				if (defined $HASH{$type}{$chrom}{$position}{2}) { $substitution{$HASH{$type}{$chrom}{$position}{2}}++; }
				print join(',',keys %substitution), "\n";
			}
		}
		close (OUTPUT);
	}
}
