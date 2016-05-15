package BamData;

use strict;
use RegularList;
use ComplementBase;
our @ISA = "RegularList";

sub new {
	my ($invocant,$filename,$minqual) = @_;
	my $class = ref($invocant) || $invocant;
	
	my $self = $class->SUPER::new($filename);
	
	return $self;
}

sub get_readname {
	my $self = shift;

	return $self->get_column(0);
}

sub get_flag {
	my $self = shift;
	
	return $self->get_column(1);
}

sub get_chrom {
	my $self = shift;
	
	return $self->get_column(2);
}

sub get_position {
	my $self = shift;
	
	return $self->get_column(3);
}

sub get_mapqual {
	my $self = shift;

	return $self->get_column(4);
}

sub get_cigar {
	my $self = shift;

	return $self->get_column(5);	
}

sub get_chrom_mate {
	my $self = shift;
	
	return $self->get_column(6);
}

sub get_position_mate {
	my $self = shift;
	
	return $self->get_column(7);
}

sub get_tlen {
	my $self = shift;
	
	return $self->get_column(8);
}

sub get_seq {
	my $self = shift;
	
	return $self->get_column(9);
}

sub get_qual {
	my $self = shift;
	
	return $self->get_column(10);
}
	
1;
