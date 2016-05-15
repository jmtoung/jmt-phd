package GetQstatNum;

require Exporter;

use strict;

our @ISA = qw(Exporter);
our @EXPORT = qw(complement_base);

sub qstat_wc_l {
	my $num = `qstat | wc -l`;
}

1;
