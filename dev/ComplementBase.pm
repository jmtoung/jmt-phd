package ComplementBase;

require Exporter;

use strict;

our @ISA = qw(Exporter);
our @EXPORT = qw(complement_base);

sub complement_base {
	my @bases = split('',shift);
	my $bases;
	
	foreach my $base (@bases) {
		if ($base eq 'A') { $bases .= 'T'; }
		elsif ($base eq 'a') { $bases .= 't'; }
		elsif ($base eq 'C') { $bases .= 'G'; }
		elsif ($base eq 'c') { $bases .= 'g'; }
		elsif ($base eq 'T') { $bases .= 'A'; }
		elsif ($base eq 't') { $bases .= 'a'; }
		elsif ($base eq 'G') { $bases .= 'C'; }
		elsif ($base eq 'g') { $bases .= 'c'; }
		else { $bases .= $base; }
	}
	return $bases if $bases;
	return undef;
}

1;
