#!/usr/bin/perl -w

use strict;

my $currentReadName;

while(<STDIN>) {
	chomp;
	
	my @split = split('\t');
	
	unless (isEqual($currentReadName,$split[0])) {
		$currentReadName = $split[0];
		
		my $seq = $split[9];
		my $qual = $split[10];
		
		my $flag = $split[1];
		
		my $strand = getStrand($flag);

		if ($strand eq '+') {
		
		} elsif ($strand eq '-') {
			$seq = reverse($seq);
			$seq =~ tr/ACGT/TGCA/;
			$qual = reverse($qual);		
		} else {
			die "[STDERR]: weird strand $strand\n";
		}
		
		my @entry = ('@' . $split[0] . '#0/1',$seq,'+',$qual);
		
		print join("\n",@entry), "\n";
	
	}
}


sub isEqual {
	my ($current,$split) = @_;
	
	return 0 if !defined $current;
	
	return 1 if $current eq $split;
	
	return 0;
}

sub getStrand {
	my $flag = shift;
	
	return '-' if ($flag & 16);
	return '+';
}
