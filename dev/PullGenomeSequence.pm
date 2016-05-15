package PullGenomeSequence;

require Exporter;

use strict;

our @ISA = qw(Exporter);
our @EXPORT = qw(pull_sequence);

sub pull_sequence {
        my $region = shift;
        my $index = shift;
        
        my $sequence = '';
        
	my $partial_sequence = `samtools faidx $index $region`;
	my @partial_sequence = split('\n',$partial_sequence);           
	shift(@partial_sequence); ### delete the first line
	$sequence .= join('',@partial_sequence);
        
        return $sequence;
}

1;
