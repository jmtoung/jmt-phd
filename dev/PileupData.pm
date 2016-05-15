package PileupData;

use strict;
use RegularList;
use ComplementBase;
our @ISA = "RegularList";

sub new {
	my ($invocant,$filename,$minqual) = @_;
	my $class = ref($invocant) || $invocant;
	
	my $self = $class->SUPER::new($filename);

	$self->{'minqual'} = $minqual;
	$self->parse_pileup if $self->has_next;
	
	return $self;
}

sub parse_next_line {
	my $self = shift;

	$self->update_next_line;
	delete $self->{'+'};
	delete $self->{'-'};
	$self->parse_pileup if $self->has_next;
}

sub get_chrom {
	my $self = shift;

	return $self->get_column(0);
}

sub get_position {
	my $self = shift;

	return $self->get_column(1);
}

sub get_ref_base {
	my $self = shift;

	return uc($self->get_column(2));
}

sub get_total_count {
	my ($self,$strands) = @_;

	### this gives counts after the minqual filter (as opposed to $self->get_column(3))
	my $count = 0;
	foreach my $strand (@{$strands}) {
		next unless exists $self->{$strand};
		foreach my $base (keys %{$self->{$strand}}) {
			$count += $self->{$strand}{$base}{'count'};
		}
	}
	
	return $count;
}

sub get_base_codes {
	my $self = shift;

	return $self->get_column(4);
}

sub get_qual_codes {
	my $self = shift;

	return $self->get_column(5);
}

sub get_bases {
	my ($self,$strands) = @_;

	my %bases;
	foreach my $strand (@{$strands}) {
		next unless exists $self->{$strand};
		foreach my $base (keys %{$self->{$strand}}) {
			$base = complement_base($base) if $strand eq '-' && scalar(@{$strands}) == 2;
			$bases{$base}++;
		}
	}
	my @bases = keys %bases;
	return \@bases;
}

sub get_base_count {
	my ($self,$strands,$bases) = @_;

	my $count = 0;
	foreach my $base (@{$bases}) {
		foreach my $strand (@{$strands}) {
			next unless exists $self->{$strand};
			$base = complement_base($base) if $strand eq '-' && scalar(@{$strands}) == 2;
			$count += $self->{$strand}{$base}{'count'} if exists $self->{$strand}{$base};
		}
	}
	return $count;
}

sub get_base_qual {
	my ($self,$strands,$bases) = @_;

	my @qual;
	foreach my $base (@{$bases}) {
		foreach my $strand (@{$strands}) {
			next unless exists $self->{$strand};
			$base = complement_base($base) if $strand eq '-' && scalar(@{$strands}) == 2;
			push(@qual,@{$self->{$strand}{$base}{'qual'}}) if exists $self->{$strand}{$base}{'qual'};
		}
	}
	return \@qual;
}

sub parse_pileup {
	my $self = shift;
	
	my @base_codes = split('',$self->get_base_codes);
	my @qual = @{convert_quality($self->get_qual_codes)};
	
	while (@base_codes || @qual) {
		my ($base_code,$base,$strand,$qual);
		$base_code = shift(@base_codes);
		
		### START OF READ
		if ($base_code eq '^') {
			shift(@base_codes); ### discard map quality
			next;
		}
		### END OF READ
		if ($base_code eq '$') {
			next;
		} 
		### DELETION
		### e.g. "-2at", which signals that the next upcoming 2 positions are deletions
		if ($base_code eq '-') {
			my $del_length;
			while($base_codes[0] =~ /[0-9]/) { $del_length .= shift(@base_codes); }
			### store in the $self->{'deletion'} hash the bases the positions, bases, and counts of the deletion
			for (my $i = 1; $i <= $del_length; $i++) {
				my $del_position = $self->get_position + $i;
				my $del_base = shift(@base_codes); ### e.g. 'a','t' in our example
				$del_base = uc($del_base);
				$self->{'deletion'}{$del_position}{$del_base}++;
			}
			next;
		}
		
		$qual = shift(@qual);

		### DELETION
		### e.g. "*" which are the actual bases that are deleted, as opposed to '-2at', see above
		if ($base_code eq '*') {
			### get one of the bases that are deleted; should equal ref_base (check later)
			my @del_base = keys %{$self->{'deletion'}{$self->get_position}};
			### decrease count
			$self->{'deletion'}{$self->get_position}{$del_base[0]}--;
			### delete if count is 0
			delete $self->{'deletion'}{$self->get_position}{$del_base[0]} if $self->{'deletion'}{$self->get_position}{$del_base[0]} == 0;
			delete $self->{'deletion'}{$self->get_position} if scalar(keys %{$self->{'deletion'}{$self->get_position}}) == 0;
			$base = 'X';
			my $del_base; ### the base that's deleted. check that it's equal to ref_base
			if ($del_base[0] =~ /([A|C|T|G|N])/) {
				$del_base = $1;
				$strand = '+';
			} elsif ($del_base[0] =~ /([a|c|t|g|n])/) {
				$del_base = uc($1);
				$strand = '-';
			}
			$del_base eq $self->get_ref_base or die "[STDERR]: error in deletion ref_base; line => ", join("\t",@{$self->get_next_line}), "\n";
		### POSITIVE STRAND MATCH 
		} elsif ($base_code eq '.') {
			$base = $self->get_ref_base;
			$strand = '+';
		### NEGATIVE STRAND MATCH
		} elsif ($base_code eq ',') {
			$base = complement_base($self->get_ref_base);
			$strand = '-';
		} elsif ($base_code =~ /([A|C|T|G|N])/) { ### POSITIVE STRAND MISMATCH
			$base = $1;
			$strand = '+';
		} elsif ($base_code =~ /([a|c|t|g|n])/) { ### NEGATIVE STRAND MISMATCH
			$base = uc(complement_base($1));
			$strand = '-';
		} elsif ($base_code eq '>') { ### SPLICE READ ON POSITIVE STRAND
			$base = 'S';
			$strand = '+';
		} elsif ($base_code eq '<') { ### SPLICE READ ON NEGATIVE STRAND
			$base = 'S';
			$strand = '-';
		}
		
		### INSERTION
		### check to see if we need to tack on anything to the base (insertion)
		### note: need to check if there is even anything inside @base_codes (might be at the end)
		if (@base_codes && $base_codes[0] eq '+') {
			shift(@base_codes); ### shift out the '+'
			my $ins_length;
			while($base_codes[0] =~ /[0-9]/) { $ins_length .= shift(@base_codes); } ### shift to get the length of the indel
			### add bases to $base
			while($ins_length-- > 0) {
				my $ins_base_code = shift(@base_codes);
				$base .= uc($ins_base_code) if $strand eq '+';
				$base .= uc(complement_base($ins_base_code)) if $strand eq '-';
			}
		}
				
		### filter by Phred quality score
		next if $self->{'minqual'} && $qual < $self->{'minqual'};

		### add pileup data to object
		$base = uc($base);
		$self->{$strand}{$base}{'count'}++;
		push(@{$self->{$strand}{$base}{'qual'}},$qual);
	}
}

sub convert_quality {
	my $temp_qual = shift;
	
	my @final_qual;
	foreach my $bit (split('',$temp_qual)) {
		#my $qual = ord($bit);
		#push(@final_qual,$qual);
		push(@final_qual,ord($bit) - 33);
	}
	
	return \@final_qual;
}

1;
