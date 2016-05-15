package PhredQual;

require Exporter;

use strict;

our @ISA = qw(Exporter);
our @EXPORT = qw(convertAscii2Phred convertPhred2Ascii);

sub convertAscii2Phred {
        my $ascii = shift;

	return ord($ascii) - 33;
}

sub convertPhred2Ascii {
	my $phred = shift;
	
	return chr($phred + 33);
}
                                               	
1;
