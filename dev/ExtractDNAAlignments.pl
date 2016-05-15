#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use lib '/ifs/h/toung/dev','/home/jmtoung/Lab/dev';

my $tag;
my $uniqueOnly = 0;
my $validOnly = 1;

$|++;

my $option = GetOptions(
	"tag=s" => \$tag,
	"uniqueOnly=i" => \$uniqueOnly,
	"validOnly=i" => \$validOnly
);

defined $tag or die "[STDERR]: define tag\n";
print STDERR "tag:\t$tag\n";
($uniqueOnly && !$validOnly) || (!$uniqueOnly && $validOnly) or die "[STDERR]: unique and valid only parameters set improperly\n";
print STDERR "uniqueOnly:\t$uniqueOnly\n";
print STDERR "validOnly:\t$validOnly\n"; 

my %current;
my %stats;
while(<STDIN>) {
	chomp;
	
	if (/^@/) {
		print $_, "\n";
		next;
	}

	my @split = split('\t');
	
	unless (is_equal(\@split,\%current)) {
		print_stuff(\%current,\%stats,$uniqueOnly,$validOnly);
		%current = ();
	}

	add_to_current(\@split,\%current,\%stats);
}

print_stuff(\%current,\%stats,$uniqueOnly,$validOnly);


print STDERR "UploadGsnapPairedBam2DB_stats", "\t", $tag, "\t", "total_num_valid_aln", "\t", $stats{'valid_aln'}{'total'}, "\n";
foreach my $freq (sort {$a<=>$b} keys %{$stats{'valid_aln'}{'freq'}}) {
	print STDERR "UploadGsnapPairedBam2DB_stats", "\t", $tag, "\t", "valid_aln_freq", "\t", $freq, "\t", $stats{'valid_aln'}{'freq'}{$freq}, "\n";
}

foreach my $flag (sort {$a<=>$b} keys %{$stats{'tally_flags'}}) {
	print STDERR "UploadGsnapPairedBam2DB_stats", "\t", $tag, "\t", "tally_flags", "\t", $flag, "\t", $stats{'tally_flags'}{$flag}, "\n";
}

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
	my ($current,$stats,$uniqueOnly,$validOnly) = @_;
	
	my @flags_counts; ### array with each flag and the counts
	foreach my $flag (sort {$a<=>$b} keys %{$current->{'flags'}}) {
		my $key = $flag . "-" . $current->{'flags'}{$flag};
		push(@flags_counts,$key);
	}
	
	my $flags_counts = join(",",@flags_counts); ### this is is so it's in order of flag
	
	my $uniqueAln = 0;
	my $validAln = 0;
	if ((defined $current->{'flags'}{99} && defined $current->{'flags'}{147}) || (defined $current->{'flags'}{83} && defined $current->{'flags'}{163})) {
		$validAln = 1;
		
		$stats->{'valid_aln'}{'total'}++;

		### this means one unique alignment
		if ($flags_counts eq '99-1,147-1' || $flags_counts eq '83-1,163-1') {
			$uniqueAln = 1;
			$stats->{'valid_aln'}{'freq'}{1}++;
		} else {
			my $secondary_count = 0;
			
			foreach my $flags ([355,403],[339,419]) {

				my $flag_count = get_flags_count($current->{'flags'},$flags);
				
				### if length of @{$flag_count} is greater than 1, that just means there's different number of flags counts, and we just take the min which is $flag_count->[0]!
				$secondary_count += $flag_count->[0] if scalar(@{$flag_count}) >= 1;
			}
			
			$stats->{'valid_aln'}{'freq'}{$secondary_count}++;
		}
		
	}

	if ($uniqueOnly) {
		return if !$uniqueAln;
	} elsif ($validOnly) {
		return if !$validAln;
	}

	foreach my $aln (@{$current->{'aln'}}) {
		print join("\t",@{$aln}), "\n";
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

