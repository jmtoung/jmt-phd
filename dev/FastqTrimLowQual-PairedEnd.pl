#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

my $fastq1;
my $fastq2;
my $tag;
my $min_length;
my $removeFailed = 1; ### added on 3/27/2012
my $clusterTag = 1;

my $result = GetOptions (
	"fastq1=s" => \$fastq1,
	"fastq2=s" => \$fastq2,
	"tag=s" => \$tag,
	"min_length=i" => \$min_length,
	"removeFailed=i" => \$removeFailed,
	"clusterTag=i" => \$clusterTag
);

my ($fastq1_name, $fastq1_dir, $fastq1_ext) = fileparse($fastq1,'\.fastq(\.gz)?');
my ($fastq2_name, $fastq2_dir, $fastq2_ext) = fileparse($fastq2,'\.fastq(\.gz)?');

$min_length =~ /[0-9]+/ or die "[STDERR]: invalid min length\n";

### check fastq is absolute name
substr($fastq1_dir,0,1) eq '/' && substr($fastq2_dir,0,1) or die "[STDERR]: '$fastq1' or '$fastq2' is not absolute file name\n";

my ($fastq1_fh, $fastq2_fh);
if ($fastq1_ext eq '.fastq') { $fastq1_fh = $fastq1; }
elsif ($fastq1_ext eq '.fastq.gz') { $fastq1_fh = "gunzip -c $fastq1 |"; }
else { die "[STDERR]: $fastq1 doesn't end in '.fastq' or '.fastq.gz'"; }
open(FASTQ1,$fastq1_fh) or die "[STDERR]: cannot open '$fastq1_fh': $!\n";

if ($fastq2_ext eq '.fastq') { $fastq2_fh = $fastq2; }
elsif ($fastq2_ext eq '.fastq.gz') { $fastq2_fh = "gunzip -c $fastq2 |"; }
else { die "[STDERR]: $fastq2 doesn't end in '.fastq' or '.fastq.gz'"; }
open(FASTQ2,$fastq2_fh) or die "[STDERR]: cannot open '$fastq2_fh': $!\n";

my $fastq1_new = $fastq1_dir . $fastq1_name . "_" . $tag . $fastq1_ext;
$fastq1_new =~ s/\.gz//;
open(FASTQ1_NEW,'>'.$fastq1_new) or die "[STDERR]: cannot open '$fastq1_new': $!\n";

my $fastq2_new = $fastq2_dir . $fastq2_name . "_" . $tag . $fastq2_ext;
$fastq2_new =~ s/\.gz//;
open(FASTQ2_NEW,'>'.$fastq2_new) or die "[STDERR]: cannot open '$fastq2_new': $!\n";

my %stats;
while(my $line1_1 = <FASTQ1>) {
	defined (my $line2_1 = <FASTQ1>) or die "[STDERR]: fastq1 line2 (sequence) undefined\n";
	defined (my $line3_1 = <FASTQ1>) or die "[STDERR]: fastq1 line3 (read_name) undefined\n";
	defined (my $line4_1 = <FASTQ1>) or die "[STDERR]: fastq1 line4 (qual scores) undefined\n";

	defined (my $line1_2 = <FASTQ2>) or die "[STDERR]: fastq2 line1 (read_name) undefined\n";
	defined (my $line2_2 = <FASTQ2>) or die "[STDERR]: fastq2 line2 (sequence) undefined\n";
	defined (my $line3_2 = <FASTQ2>) or die "[STDERR]: fastq2 line3 (read_name) undefined\n";
	defined (my $line4_2 = <FASTQ2>) or die "[STDERR]: fastq2 line4 (qual scores) undefined\n";

	$stats{'number_total_sequenced'}++;

	chomp ($line1_1, $line2_1, $line3_1, $line4_1, $line1_2, $line2_2, $line3_2, $line4_2);

	my @line1_1 = split('\s+',$line1_1);
	my @line1_2 = split('\s+',$line1_2);

	### this part added on 3/27/2012
	my @line1_info_1 = split(':',$line1_1[1]);
	my @line1_info_2 = split(':',$line1_2[1]);
	
	if ($clusterTag && $removeFailed) {
		($line1_info_1[1] eq 'N' or $line1_info_1[1] eq 'Y') or die "[STDERR]: cluster tag doesn't equal Y/N: $line1_1\n";
		($line1_info_2[1] eq 'N' or $line1_info_2[1] eq 'Y') or die "[STDERR]: cluster tag doesn't equal Y/N: $line1_2\n";
		 
		$stats{'removeFailed'}{$line1_info_1[1]}{$line1_info_2[1]}++;

		next unless $line1_info_1[1] eq 'N' && $line1_info_2[1] eq 'N';
	}
	###
	
	my @line2_1 = split('',$line2_1);
	my @line2_2 = split('',$line2_2);
	
	my @line4_1 = split('',$line4_1);
	my @line4_2 = split('',$line4_2);

	### while @line4 is defined AND the last qual score is #, pop it off
	### @line4 might not be defined b/c the entire qual score might all be #'s
	### keep track of number popped off
	my $num_pop_1 = 0;
	my $num_pop_2 = 0;
	while(@line4_1 && $line4_1[-1] eq '#') {
		pop(@line4_1);
		pop(@line2_1);
		$num_pop_1++;
	}

	while(@line4_2 && $line4_2[-1] eq '#') {
		pop(@line4_2);
		pop(@line2_2);
		$num_pop_2++;
	}

	$line1_1[0] .= "|" . $num_pop_1 . ":" . $num_pop_2; ### fixed this on 3/31/2012
	$line1_2[0] .= "|" . $num_pop_1 . ":" . $num_pop_2;
	
	$stats{'number_bases_trimmed'}{$num_pop_1}{$num_pop_2}++;

	### if read (@line2) length and quality score (@line4) length are not the same, exit
	if ((scalar(@line2_1) != scalar(@line4_1)) || (scalar(@line2_2) != scalar(@line4_2))) {
		die "[STDERR]: 'sequence' and 'quality score' length unequal =>\n$line1_1\n$line2_1\n$line3_1\n$line4_1\n$line1_2\n$line2_2\n$line3_2\n$line4_2\n";
	} 

	### check if read lengths >= min_length	
	my $fastq1_ge_min = 0;
	my $fastq2_ge_min = 0;
	$fastq1_ge_min = 1 if (scalar(@line2_1) >= $min_length);
	$fastq2_ge_min = 1 if (scalar(@line2_2) >= $min_length);
	

	### if the length is not greater than min length, then skip
	$stats{'ge_min_length'}{$fastq1_ge_min}{$fastq2_ge_min}++;
	next unless $fastq1_ge_min && $fastq2_ge_min;

	print FASTQ1_NEW join(" ",@line1_1), "\n", join('',@line2_1), "\n", $line3_1, "\n", join('',@line4_1), "\n";
	print FASTQ2_NEW join(" ",@line1_2), "\n", join('',@line2_2), "\n", $line3_2, "\n", join('',@line4_2), "\n";
}
close(FASTQ1);

my $see_if_remain = <FASTQ2>;
if (defined $see_if_remain) {
	die "[STDERR]: fastq2 still has remaining lines!\n";
}

my @report = ($fastq1,$fastq2,$fastq1_new,$fastq2_new);
print STDOUT "stats", "\t", "number_total_sequenced", "\t", join("\t",@report), "\t", $stats{'number_total_sequenced'}, "\n";

foreach my $num1 (sort {$a<=>$b} keys %{$stats{'number_bases_trimmed'}}) {
	foreach my $num2 (sort {$a<=>$b} keys %{$stats{'number_bases_trimmed'}{$num1}}) {
		print STDOUT "stats", "\t", "number_bases_trimmed", "\t", join("\t",@report), "\t", $num1, "\t", $num2, "\t", $stats{'number_bases_trimmed'}{$num1}{$num2}, "\n";
	}
}
foreach my $ge1 (sort {$a<=>$b} keys %{$stats{'ge_min_length'}}) {
	foreach my $ge2 (sort {$a<=>$b} keys %{$stats{'ge_min_length'}{$ge1}}) {
		print STDOUT "stats", "\t", "number_ge_min_length", "\t", join("\t",@report), "\t", $ge1, "\t", $ge2, "\t", $stats{'ge_min_length'}{$ge1}{$ge2}, "\n";
	}
}

if ($removeFailed) {
	foreach my $fail1 (sort keys %{$stats{'removeFailed'}}) {
		foreach my $fail2 ( sort keys %{$stats{'removeFailed'}{$fail1}}) {
			print STDOUT "stats", "\t", "removeFailed", "\t", join("\t",@report), "\t", $fail1, "\t", $fail2, "\t", $stats{'removeFailed'}{$fail1}{$fail2}, "\n";
		}
	}
}

print "completed_successfully\t", join("\t",@report), "\n";

