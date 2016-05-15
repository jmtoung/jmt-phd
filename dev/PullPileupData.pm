package PullPileupData;

use strict;
use lib '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev';
use PileupData;
use ComplementBase;

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = {};
	bless $self, $class;

	$self->{'bam'} = shift;
	$self->{'index'} = shift;
	$self->{'minqual'} = shift;
	-e $self->{'bam'} && -e $self->{'index'} && defined $self->{'minqual'} or die "[STDERR]: undefined values in PullPileupData\n";	
	return $self;
}

sub get_pileup {
	my $self = shift;
	my $region = shift;
	
	my $bam = $self->{'bam'};
	my $index = $self->{'index'};
	my $minqual = $self->{'minqual'};

	my $pileup = "samtools mpileup -f $index -r $region $bam |";
	open(PILEUP,$pileup) or die "[STDERR]: can't fork $pileup: $!\n";
	my @results;
	while(<PILEUP>) {
		chomp;
		push(@results,$_);
	}
	
	return PileupData->new(\@results,$minqual);
}

1;
