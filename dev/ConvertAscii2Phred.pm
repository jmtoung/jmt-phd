package ConvertAscii2Phred;

require Exporter;

use strict;

our @ISA = qw(Exporter);
our @EXPORT = qw(convert_quality);

sub convert_quality {
        my $temp_qual = shift;
        
	my @final_qual;
	foreach my $bit (split('',$temp_qual)) {
		push(@final_qual,ord($bit) - 33);
	}                                                
	return \@final_qual;
}
                                                	
1;
