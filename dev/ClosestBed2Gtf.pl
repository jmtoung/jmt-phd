#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Tie::IxHash;
use File::Basename;
use lib '/home/jmtoung/Lab/dev';
use List::Util qw[min max];

my $bed;
my $gtf_files;

my $options = GetOptions(
	"bed=s" => \$bed,
	"gtf_files=s" => \$gtf_files
);

my @gtf_files = split(',',$gtf_files);
foreach my $gtf_file (@gtf_files) {
	-e $gtf_file or die "[STDERR]: gtf file $gtf_file doesn't exist\n";
}
$gtf_files = join(' ',@gtf_files);

my $closestBed = "cat $gtf_files | closestBed -d -a $bed -b stdin |";

open(COMMAND,$closestBed) or die "[STDERR]: can't open $closestBed: $!\n";

my %current; ### the current entry we're working on
while(<COMMAND>) {
	chomp;
	my @split = split('\t'); ### a line of output from 'closestBed' command
	
	if (is_equal(\@split,\%current)) {
		add_to_current(\@split,\%current);
	} else {	
		calculate_stuff(\%current);
		%current = ();
		add_to_current(\@split,\%current);
	}
	
}
calculate_stuff(\%current);
close(COMMAND);

sub is_equal {
	my ($split,$current) = @_;
	
	return 1 if !defined $current->{'chrom'};
	### equality means chromosome, start, end, and strand are all equal. (the input to closestBed)
	return 0 unless $current->{'chrom'} eq $split->[0];
	return 0 unless $current->{'start'} == $split->[1];
	return 0 unless $current->{'end'} == $split->[2];
	return 0 unless $current->{'strand'} eq $split->[3];
	
	return 1;
}

sub add_to_current {
	my ($split,$current) = @_;

	### define chromosome, start, end, strand; note that you can overwrite if they are already defined since we checked they're equal
	$current->{'chrom'} = shift(@{$split});
	$current->{'start'} = shift(@{$split});
	$current->{'end'} = shift(@{$split});
	$current->{'strand'} = shift(@{$split});

	### make a hash that contains all other information	
	my %gene_record;
	$gene_record{'type'} = $split->[2];
	$gene_record{'start'} = $split->[3];
	$gene_record{'end'} = $split->[4];
	$gene_record{'strand'} = $split->[6];
	foreach my $attribute (split(';',$split->[8])) { if ($attribute =~ /(\w+)\s+"(.*)"/) { $gene_record{$1} = $2; } }
	$gene_record{'dist'} = $split->[9];

	push(@{$current->{'gene_records'}},\%gene_record);
	
}

sub calculate_stuff {
	my ($current,$closestBed) = @_;

	my @line = ($current->{'chrom'},$current->{'start'},$current->{'end'},$current->{'strand'});
	my %dist; ### dist of entry (A) to closest entry in B
	my %feature_types; ### e.g. exon, intron, etc.
	my %gene_info; ### stuff like transcript_id, etc.
	my @info_types = ('transcript_id','gene_name','gene_type','transcript_type');
	### calculate type of RDD by looping through all the "gene_record" hashes
	foreach my $gene_record (@{$current->{'gene_records'}}) {
		$dist{$gene_record->{'dist'}}++;
		### store info for this gene
		my $gene_id = $gene_record->{'gene_id'};
		foreach my $info_type (@info_types) { 
			my $info = $gene_record->{$info_type};
			next if $info_type eq 'transcript_id' && $info =~ /ENSG/;
			$gene_info{$gene_id}{$info_type}{$info}++ if defined $info;

		}
		
		### see if rdd falls within record bounds to determine the type; if it doesn't fall w/in, it's intergenic
		if ($current->{'end'} <= $gene_record->{'end'} && $current->{'end'} >= $gene_record->{'start'}) {
			### suppress 'gene' or 'transcript' types b/c they're redundant info...
			$feature_types{$gene_record->{'type'}}++ unless $gene_record->{'type'} eq 'gene' or $gene_record->{'type'} eq 'transcript';
		} else {
			$feature_types{'intergenic'}++;
		}
	}
	
	### calculate_gene_info
	### first aggregate gene info by gene_id
	my %agg_gene_info;
	tie %agg_gene_info, "Tie::IxHash";
	foreach my $gene_id (sort keys %gene_info) {
		$agg_gene_info{'gene_id'}{$gene_id}++ unless $gene_id =~ /^(NM|NR)_[0-9]+/;
		foreach my $info_type (@info_types) {
			my $info = join(",",sort keys %{$gene_info{$gene_id}{$info_type}});
			$agg_gene_info{$info_type}{$info}++ unless $info eq '';
		}
	}

	### push the aggregated gene info across all genes onto the @line
	foreach my $info_type ('gene_id',@info_types) {
		my $INFO = join(";",keys %{$agg_gene_info{$info_type}});
		$INFO = '-' if $INFO eq '';
		push(@line,$INFO);
	}

	### if no gene record fell directly on this RDD site, this means it's intronic.
	### if a gene record fell directly on this RDD site, it should be exonic
	### if nothing fell, and the closest thing did not envelop RDD, then it's intergenic
	if (scalar(keys %feature_types) == 0) { $feature_types{'intron'}++; }
	push(@line,join(",",sort keys %feature_types));
	
	### calculate distance from transcription start site
	my (@gene_start,@gene_end,%gene_strand);
	foreach my $gene_record (@{$current->{'gene_records'}}) {
		if ($gene_record->{'type'} eq 'gene') {
			push(@gene_start,$gene_record->{'start'});
			push(@gene_end,$gene_record->{'end'});
		}
		$gene_strand{$gene_record->{'strand'}}++;
	}
	my $gene_strand = join(",",sort keys %gene_strand);	
	
	my $distance;
	if (@gene_start == 1 && $gene_strand ne '+,-') {
		if ($gene_strand eq '-') { $distance = $gene_end[0] - $current->{'end'}; }
		elsif ($gene_strand eq '+') { $distance = $current->{'end'} - $gene_start[0]; }
		else { die "[STDERR]: gene_strand weird, can't define distance\n"; }
		push(@line,$distance);
	} else {
		push(@line,'NA');
	}

	push(@line,$gene_strand);

	my $dist = min(keys %dist);
	push(@line,$dist);
	
	print join("\t",@line), "\n";
}
