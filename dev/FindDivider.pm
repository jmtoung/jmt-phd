package FindDivider;

require Exporter;

use strict;

our @ISA = qw(Exporter);
our @EXPORT = qw(find_divider);

sub find_divider {
	my $split = shift;
	
	for (my $i = 0; $i < scalar(@{$split}); $i++) {
		if ($split->[$i] =~ /^total_count;/) {
			return $i;
		}
	}
	return undef;
}

1;
