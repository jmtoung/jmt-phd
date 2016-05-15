package OrderedGenomicList;

use IO::File;
use strict;
use File::Basename;

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = { };
	bless $self, $class;
	
	my $filename = shift;
	open($self->{'file'},$filename) or die "couldn't open $filename: $!\n";
	
	($self->{'filename'},$self->{'filedir'},$self->{'fileext'}) = fileparse($filename,'\.[a-zA-Z0-9]+$');

	if(shift) { $self->{'next_line'} = $self->read_line; }
	else { 	$self->{'next_line'} = ['placeholder']; }### this is the placeholder before first line is returned (just to say that we do indeed have a next line, which is true since open(file) was successful)
	
	return $self;
}

### get methods

sub get_file {
	my $self = shift;
	return $self->{'file'};
}

sub get_filename {
	my $self = shift;
	return $self->{'filename'};
}

sub get_filedir {
	my $self = shift;
	return $self->{'filedir'};
}

sub get_fileext {
	my $self = shift;
	return $self->{'fileext'};
}

sub get_next_line {
	my $self = shift;
	return $self->{'next_line'};
}

### other methods

sub update_next_line {
	my $self = shift;
	$self->{'next_line'} = $self->read_line();
}

sub read_line {
	my $self = shift;
	my $fh = $self->get_file();
	my $line = <$fh>;
	if (defined $line) { 
		chomp $line;
		my @line = split('\t',$line); 
		return \@line;
	}
	return $line;
}

sub has_next {
	my $self = shift;
	return defined $self->get_next_line;
}

sub get_column {
	my $self = shift;
	my $column = shift;
	if ($self->has_next()) {
		return $self->get_next_line->[$column];
	} else {
		return $self->get_next_line;
	}
}

1;
