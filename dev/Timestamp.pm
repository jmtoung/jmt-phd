package Timestamp;

require Exporter;

use strict;

our @ISA = qw(Exporter);
our @EXPORT = qw(get_timestamp);

sub get_timestamp {
	my $TIME = localtime;
	my @TIME = split('\s+',$TIME);
	
	### get rid of day of week
	shift(@TIME);

	### convert from Jan to January
	$TIME[0] = full_month($TIME[0]);
	
	return join(" ",@TIME);
}

sub full_month {
	my $MONTH = shift;
	return 'January' if $MONTH eq 'Jan';
	return 'February' if $MONTH eq 'Feb';
	return 'March' if $MONTH eq 'Mar';
	return 'April' if $MONTH eq 'Apr';
	return 'May' if $MONTH eq 'May';
	return 'June' if $MONTH eq 'Jun';
	return 'July' if $MONTH eq 'Jul';
	return 'August' if $MONTH eq 'Aug';
	return 'September' if $MONTH eq 'Sep';
	return 'October' if $MONTH eq 'Oct';
	return 'November' if $MONTH eq 'Nov';
	return 'December' if $MONTH eq 'Dec';
	return $MONTH;
}

1;
