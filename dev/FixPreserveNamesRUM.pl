#!/usr/bin/perl -w

use lib '/gpfs/fs121/h/toung/dev';
use RegularList;
use strict;
use Getopt::Long;

my $rumFolder;

my $options = GetOptions(
	"rumFolder=s" => \$rumFolder
);

my $map = $rumFolder . "/" . "read_names_mapping";

my $sam = $rumFolder . "/" . "RUM.sam";

-e $map && -e $sam or die "[STDERR]: $map or $sam doesn't exist: $!\n";

print STDERR "read_names_mapping:\t$map\n";
print STDERR "sam:\t$sam\n";

my $MAP = RegularList->new($map);
open(SAM,$sam) or die "[STDERR]: can't open $sam: $!\n";

while(<SAM>) {
	chomp;
	
	my @split = split('\t');
	
	if ($split[0] =~ /^@/) {
		print $_, "\n";
		next;
	}

	my $readName = $split[0];
	my @readName = split('\.',$readName);
	my $readId = $readName[1];
	
	while (!defined @{$MAP->get_next_line}) { $MAP->update_next_line; }

	while($MAP->has_next && $readId != getSeqId($MAP->get_column(0))) {

		$MAP->update_next_line;
		
		while (!defined @{$MAP->get_next_line}) { $MAP->update_next_line; }
	
		$readId >= getSeqId($MAP->get_column(0)) or die "[STDERR]: readId is greater than seqId\n";

	}
	
	$readId == getSeqId($MAP->get_column(0)) or die "[STDERR]: readId not equal to seqId\n";
	
	$split[0] = formatReadName($MAP->get_column(1));
	
	print join("\t",@split), "\n";

}
close(SAM);

sub getSeqId {
	my $readName = shift;
	
	$readName =~ s/(a|b)$//;
	
	my @readName = split('\.',$readName);
	
	return $readName[1];
}

sub formatReadName {
	my $readName = shift;
	
	$readName =~ s/\/(1|2)$//;
	
	return $readName;
	
}
