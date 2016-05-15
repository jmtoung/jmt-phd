package IlmnDNA;

use strict;
use RegularList;
use ComplementBase;
our @ISA = "RegularList";

sub new {
	my ($invocant,$filename) = @_;
	my $class = ref($invocant) || $invocant;
	
	my $self = $class->SUPER::new($filename);

	return $self;
}

sub get_chrom_NCBI {
	my $self = shift;
	
	my $chrom = $self->get_column(0);
	$chrom =~ s/chr//;
	$chrom = 'MT' if $chrom eq 'M';
	
	return $chrom;
}

sub get_start {
	my $self = shift;
	
	return $self->get_column(1);
}

sub get_end {
	my $self = shift;
	
	return $self->get_column(2);
}

sub get_var_type {
	my $self = shift;
	
	return $self->get_column(3);
}

sub get_depth {
	my $self = shift;
	
	return $self->get_column(4);
}

sub get_max_gt {
	my $self = shift;
	
	return $self->get_column(5);
}

sub get_max_gt_poly {
	my $self = shift;
	
	return $self->get_column(6);
}

sub get_max_gt_call {
	my $self = shift;
	
	return $self->get_column(7);
}

sub get_max_gt_poly_call {
	my $self = shift;
	
	return $self->get_column(8);
}

sub get_num_A {
	my $self = shift;
	
	return $self->get_column(9);
}

sub get_num_C {
	my $self = shift;
	
	return $self->get_column(10);
}

sub get_num_G {
	my $self = shift;
	
	return $self->get_column(11);
}

sub get_num_T {
	my $self = shift;
	
	return $self->get_column(12);
}
