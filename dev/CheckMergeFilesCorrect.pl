#!/usr/bin/perl -w

use strict;

my %data;
while(<STDIN>) {
	chomp;
	
	my @split = split('\t');
	
	my @fileName = split('\.',$split[1]);
	$fileName[-2] = $fileName[-1];
	pop(@fileName);

	my $fileName = join(".",@fileName);
	
	$data{$fileName}{$split[2]}++;	
	
}


foreach my $file (sort keys %data) {
	my $numKeys = keys %{$data{$file}};
	
	my $totalLines;
	
	foreach my $num (sort keys %{$data{$file}}) {
		$totalLines += $num * $data{$file}{$num};
	}
	
	print join("\t",'numKeys',$file,$numKeys), "\n";
	print join("\t",'totalLines',$file,$totalLines/1000000), "\n";

}
