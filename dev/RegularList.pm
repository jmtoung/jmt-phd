package RegularList;

use strict;
use File::Basename;

sub new {
	my ($invocant,$filename,$delimeter) = @_;
	my $class = ref($invocant) || $invocant;
	my $self = {};
	bless $self, $class;

	$filename =~ s/\|/\\|/g;

	if ($filename =~ /\.gz$/) {
		open($self->{'fh'},"gunzip -c $filename |") or die "[STDERR]: can't fork $filename: $!\n";
	} elsif ($filename =~ /\.bam$/) {
		open($self->{'fh'},"samtools view -h $filename |") or die "[STDERR]: can't form $filename: $!\n";
	} else {
		open($self->{'fh'},$filename) or die "[STDERR]: couldn't open $filename: $!\n";
	}

	($self->{'file_name'},$self->{'file_dir'},$self->{'file_ext'}) = fileparse($filename,'\.\w+(\.gz)?');

	if (defined $delimeter) { $self->{'delimeter'} = $delimeter; }
	else { $self->{'delimeter'} = '\t'; }

	$self->{'next_line'} = $self->update_next_line;
	
	return $self;
}

sub get_fh {
	my $self = shift;
	
	return $self->{'fh'};
}

sub get_file_name {
	my $self = shift;
	
	return $self->{'file_name'};
}

sub get_file_dir {
	my $self = shift;
	
	return $self->{'file_dir'};
}

sub get_file_ext {
	my $self = shift;
	
	return $self->{'file_ext'};
}

sub get_next_line {
	my $self = shift;
	
	return $self->{'next_line'};
}

sub update_next_line {
	my $self = shift;
	
	my $fh = $self->get_fh;
	my $line = <$fh>;
	if (defined $line) { 
		chomp $line;
		my @line = split($self->{'delimeter'},$line);
		$self->{'next_line'} = \@line;
	} else {
		$self->{'next_line'} = undef;
	}
}

sub has_next {
	my $self = shift;

	return defined $self->get_next_line;
}

sub get_column {
	my ($self,$column) = @_;

	if ($self->has_next) {
		my $next_line = $self->get_next_line;
		### check if column is with bounds of $next_line (array ref)
		if ($column >= 0 && $column < scalar(@{$next_line})) {
			return $next_line->[$column];
		} else {
			die "[STDERR]: array index out of bounds in get_column: column => $column; line => (", join(",",@{$next_line}), ")\n";
		}
	} else {
		return undef;
	}
}

1;
