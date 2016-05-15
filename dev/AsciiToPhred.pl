#!/usr/bin/perl -w

use strict;
use Getopt::Long;


my $string = 'GGEDHHHHBHGGHHGHGHDGGGBDGGGGGEGGEGEBBB>GGGDG;DDDGBGG8ECDDED?BDBB@DDBC?BDB@HHF8GHHH@DG<@>>?@?';
my $bases = 'T.T,..T,T.T.T..,TT.,,';
my $options = GetOptions(
	"string=s" => \$string
);

my $converted = convert_quality($string);
my @bases = split('',$bases);

my $index = 1;
foreach my $convert (@{$converted}) {
	print $index, "\t", $convert, "\n";
	$index++;
}

sub convert_quality {
        my $temp_qual = shift;
        
        my @final_qual;
        foreach my $bit (split('',$temp_qual)) {
                push(@final_qual,ord($bit) - 33);
        }
        
        return \@final_qual;
}



