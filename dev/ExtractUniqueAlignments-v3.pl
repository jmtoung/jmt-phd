#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use lib '/ifs/h/toung/dev','/home/jmtoung/Lab/dev';

### this script differs from the original "ExtractUniqueAlignments.pl" only in that it accepts the following flags (161,81,145,97) in addition to the previous (163,83,147,99). 
### this is b/c Tophat doesn't mark the "proper pair" bit (hence the difference of 2) sometimes accurately

### FOR PAIRED-END READS
### valid flag pairs are (163 & 83) and (147 & 99) along with the same flags - 2. (4 total)
### 'validOnly' will print a read and all its alignments if at least one valid flag pair is defined
### 'uniqueOnly' will only print a read and its alignments if the only flags defined are one of the 4 valid pairs
### 'primaryOnly' will print a read and all of the alignments that are of the 8 allowable flags if at least one valid flag pair is defined

### FOR SINGLE-END READS
### valid flags are 0, 16, which are primary alignments and 256 and 272 which are secondary alignments. 4 means nothing aligned
### 'validOnly' will print everything actually...except for 4
### 'uniqueOnly' will only print a read if its only flag is 0 or 16
### 'primaryOnly' will print only 0 and 16 for every read pair if 0 or 16 is defined (256 and 272 can be defined but wont be printed).

my $bam;
my $uniqueOnly = 0;
my $validOnly = 0;
my $primaryOnly = 1;

$|++;

my $option = GetOptions(
	"bam=s" => \$bam,
	"uniqueOnly=i" => \$uniqueOnly,
	"validOnly=i" => \$validOnly,
	"primaryOnly=i" => \$primaryOnly
);

(defined $bam && -e $bam) && print STDERR "bam:\t$bam\n" or die "[STDERR]: bam not defined\n";
my $optionsSum = $uniqueOnly + $validOnly + $primaryOnly;
$optionsSum == 0 || $optionsSum == 1 or die "[STDERR]: unique and valid only parameters set improperly\n";
print STDERR "uniqueOnly:\t$uniqueOnly\n";
print STDERR "validOnly:\t$validOnly\n"; 
print STDERR "primaryOnly:\t$primaryOnly\n";

my $command = "samtools view -h $bam |";
open(COMMAND,$command) or die "[STDERR]: can't open $command: $!\n";

my %current;
my %stats;
while(<COMMAND>) {
	chomp;
	
	if (/^@/) {
		print $_, "\n";
		next;
	}

	my @split = split('\t');
	
	unless (is_equal(\@split,\%current)) {
		print_stuff(\%current,\%stats,$uniqueOnly,$validOnly,$primaryOnly);
		%current = ();
	}

	add_to_current(\@split,\%current,\%stats);
}

print_stuff(\%current,\%stats,$uniqueOnly,$validOnly,$primaryOnly);

print STDERR "UploadPairedBam2DB_stats", "\t", $bam, "\t", "total_num_seen_reads", "\t", $stats{'seen_reads'}{'total'}, "\n";

print STDERR "UploadPairedBam2DB_stats", "\t", $bam, "\t", "total_num_valid_aln", "\t", $stats{'valid_aln'}{'total'}, "\n";
foreach my $freq (sort {$a<=>$b} keys %{$stats{'valid_aln'}{'freq'}}) {
	print STDERR "UploadPairedBam2DB_stats", "\t", $bam, "\t", "valid_aln_freq", "\t", $freq, "\t", $stats{'valid_aln'}{'freq'}{$freq}, "\n";
}

foreach my $flag (sort {$a<=>$b} keys %{$stats{'tally_flags'}}) {
	print STDERR "UploadPairedBam2DB_stats", "\t", $bam, "\t", "tally_flags", "\t", $flag, "\t", $stats{'tally_flags'}{$flag}, "\n";
}

#####################################################################

sub is_equal {
	my ($split,$current) = @_;
	
	return 1 if !defined $current->{'readname'};
	return 0 unless $split->[0] eq $current->{'readname'};
	
	return 1;
}

sub add_to_current {
	my ($split,$current,$stats) = @_;

	$current->{'readname'} = $split->[0];
		
	$current->{'flags'}{$split->[1]}++;

	$stats->{'tally_flags'}{$split->[1]}++;
	
	push(@{$current->{'aln'}},$split);
}

sub print_stuff {
	my ($current,$stats,$uniqueOnly,$validOnly,$primaryOnly) = @_;
	
	$stats->{'seen_reads'}{'total'}++;
	
	my @flags_counts; ### array with each flag and the counts
	foreach my $flag (sort {$a<=>$b} keys %{$current->{'flags'}}) {
		my $key = $flag . "-" . $current->{'flags'}{$flag};
		push(@flags_counts,$key);
	}
	
	my $flags_counts = join(",",@flags_counts); ### this is is so it's in order of flag

	my $uniqueAln = 0;
	my $validAln = 0;
	if ((defined $current->{'flags'}{99} && defined $current->{'flags'}{147}) || 
		(defined $current->{'flags'}{83} && defined $current->{'flags'}{163}) || 
		(defined $current->{'flags'}{97} && defined $current->{'flags'}{145}) || 
		(defined $current->{'flags'}{81} && defined $current->{'flags'}{161}) ||
		(defined $current->{'flags'}{0}) ||
		(defined $current->{'flags'}{16})) {

		### for PE: valid means one of the 4 flag pairs exists in combo (99 and 147, or 83 and 163, etc)
		### for SE: valid means either 0 or 16 exist	
		$validAln = 1; 
		$stats->{'valid_aln'}{'total'}++;

		### this means one unique alignment
		if ($flags_counts eq '99-1,147-1' || 
			$flags_counts eq '83-1,163-1' || 
			$flags_counts eq '97-1,145-1' || 
			$flags_counts eq '81-1,161-1' ||
			$flags_counts eq '0-1' ||
			$flags_counts eq '16-1') {

			$uniqueAln = 1;
			$stats->{'valid_aln'}{'freq'}{1}++;
			
		} else {
			my $secondary_count = 0;
			
			foreach my $flags ([355,403],[339,419],[353,401],[337,417],[256],[272]) {

				my $flag_count = get_flags_count($current->{'flags'},$flags);
				
				### if length of @{$flag_count} is greater than 1, that just means there's different number of flags counts, and we just take the min which is $flag_count->[0]!
				$secondary_count += $flag_count->[0] if scalar(@{$flag_count}) >= 1;
				
			}

			$stats->{'valid_aln'}{'freq'}{$secondary_count}++;
		}
		
	}

	if ($primaryOnly) {
	
		return unless $validAln;
		
		foreach my $aln (@{$current->{'aln'}}) {
			my $flag = $aln->[1];
			
			foreach my $ok (99,147,83,163,0,16) {
				if ($flag == $ok || ($flag == ($ok - 2))) {
					print join("\t",@{$aln}), "\n";
				}
			}
		}
	
	} else {
		if ($uniqueOnly) {
			return if !$uniqueAln;
		} elsif ($validOnly) {
			return if !$validAln;
		}
				
		foreach my $aln (@{$current->{'aln'}}) {
			print join("\t",@{$aln}), "\n";
		}	
	}
}




### this gives you the number of valid alignments (primary and then not primary)

sub calc_num_aln {
	my $hash = shift;
	
	my $primary_aln = 0;
	foreach my $flags ([99,147],[83,163]) {
		my $flags_count = flags_count($hash,$flags);
		### if both flags are defined and their counts are the same, then...
		if (flags_defined($hash,$flags) && scalar(@{$flags_count}) == 1) {
			$primary_aln += $flags_count->[0];
		}
	}
	
	my $secondary_aln = 0;
	foreach my $flags ([355,403],[339,419]) {
		my $flags_count = flags_count($hash,$flags);
		if (flags_defined($hash,$flags) && scalar(@{$flags_count}) == 1) {
			$secondary_aln += $flags_count->[0];
		}
	}

}

sub get_flags_count {
	my ($hash,$flags) = @_;
	
	my %count;
	
	foreach my $flag (@{$flags}) {
		$count{$hash->{$flag}}++ if defined $hash->{$flag};
	}
	
	my @count = sort {$a<=>$b} keys %count;
	return \@count;
}

sub flags_defined {
	my ($hash,$flags) = @_;
	
	foreach my $flag (@{$flags}) {
		return 0 if !defined $hash->{$flag};
	}
	return 1;
}


### unsplit options ############################################################
#my @options = split('\.',$options);
#my $options_split = join(" ",@options);
################################################################################

#my ($bam_name,$bam_dir,$bam_ext) = fileparse($bam,'\.bam');
#my $fastq = $bam_dir;
#$fastq =~ s/aln/fastq/;
#chop $fastq;
#$fastq = $fastq . ".fastq";
#print "fastq:\t", $fastq, "\n";

### remove $HOME from $fastq ### added 12/18/2011
#my $HOME_REGEX = $HOME;
#$HOME_REGEX =~ s/\//\\\//g;
#my $fastq_rmhome = $fastq;
#$fastq_rmhome =~ s/$HOME_REGEX//g;
#$fastq_rmhome =~ s/^\///;

### look up $fastq in fastqDB to get fastqID ###################################
#my $fastqID = Database->new($fastqDB)->lookup(1,$fastq_rmhome,0);
#defined $fastqID or die "[STDERR]: FASTQFAIL '$fastq' doesn't exist in $fastqDB\n";
################################################################################

### calculate number aligned ###################################################
#my ($num_aligned,$hist_alignment) = calculate_bam_stats($bam);
################################################################################

#my $num_warning_exhausted = 0;
#my $runtime = '-';

### remove $HOME from $bam
#my $bam_rmhome = $bam;
#$bam_rmhome =~ s/$HOME_REGEX//g;
#$bam_rmhome =~ s/^\///;

### add entry to database ######################################################
#my @entry = (	
#	$bam_rmhome,
#	$fastqID,
#	$index,
#	"gsnap",
#	sort_options($options),
#	$num_aligned,
#	$hist_alignment,
#	$num_warning_exhausted,
#	$runtime,
#	"'" . get_timestamp . "'"
#);
#print "***Adding Entry to Database ***\n";
#print "\n", join("\t",@entry), "\n";
#Database->new($alnDB)->add2db(\@entry,1);
################################################################################

sub sort_options {
	my $options = shift;

	my @options = split('\s+',$options);
	my %options;
	my $value;
	while (@options) {
		my $token = pop(@options);
		if ($token =~ /-/) { 
			$options{$token} = $value; 
			$value = undef;
		} else { 
			$value = $token; 
		}
	}

	foreach my $key (keys %options) { 
		if (defined $options{$key}) {		
			push(@options,$key . " " . $options{$key});
		} else {
			push(@options,$key);
		}
	}

	return join(" ",@options);
}

