#!/usr/bin/perl -w

use strict;
use File::Basename;

my $trash = "/home/jmtoung/Trash";

foreach my $file (@ARGV) {
	my $destination = $trash;

	$file =~ s/(\/)_$//;
	$file =~ s/\(/\\\(/g;
        $file =~ s/\)/\\\)/g;

	my ($fileName,$fileDir,$fileExt) = fileparse($file,'');
	$destination .= "/$fileName";

	my $count = 1;
	my $destinationCopy = $destination;
	while (-e $destinationCopy) {
		$destinationCopy = $destination . ".$count";
		$count++;
	}

	my $command = "mv $file $destinationCopy";

	!system($command) or die "[STDERR]: can't $command: $!\n";
}



