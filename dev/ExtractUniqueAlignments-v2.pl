#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use lib '/ifs/h/toung/dev','/home/jmtoung/Lab/dev';

#### this script takes Greg Grant's advice and simply removes any read pair where there's more than one alignment per read.

my $bam;

$|++;

my $option = GetOptions(
	"bam=s" => \$bam,
);

(defined $bam && -e $bam) && print STDERR "bam:\t$bam\n" or die "[STDERR]: bam not defined\n";

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
		print_stuff(\%current,\%stats);
		%current = ();
	}

	add_to_current(\@split,\%current,\%stats);
}

print_stuff(\%current,\%stats);

print STDERR "UploadPairedBam2DB_v2_stats", "\t", $bam, "\t", "total_num_seen_reads", "\t", $stats{'seen_reads'}{'total'}, "\n";

print STDERR "UploadPairedBam2DB_v2_stats", "\t", $bam, "\t", "total_num_unique_aln", "\t", $stats{'valid_aln'}{'total'}, "\n";

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

	my $strand = get_strand($split->[1]);
	
	$current->{$strand}++;
			
	$stats->{'tally_flags'}{$split->[1]}++;
	
	push(@{$current->{'aln'}},$split);
}

sub print_stuff {
	my ($current,$stats,$uniqueOnly,$validOnly) = @_;
	
	$stats->{'seen_reads'}{'total'}++;

	if (defined $current->{'+'} && defined $current->{'-'} && $current->{'+'} == 1 && $current->{'-'} == 1) {
	
		$stats->{'valid_aln'}{'total'}++;
	
		foreach my $aln (@{$current->{'aln'}}) {
			print join("\t",@{$aln}), "\n";
		}	
	
	}

}

sub get_strand {
	my $FLAG = shift;
	
	return '-' if ($FLAG & 16);
	return '+';
}
