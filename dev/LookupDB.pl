#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib '/ifs/h/toung/dev';
use Database;

my $database;
my $col_value;
my $value;
my $col_lookup;

my $options = GetOptions(	"database=s" => \$database,
				"col_value=i" => \$col_value, 
				"value=s" => \$value,
				"col_lookup=i" => \$col_lookup
			);

my $DB = Database->new($database);
my $value_lookup = $DB->lookup($col_value,$value,$col_lookup);

if (defined $value_lookup) {
	print $value_lookup, "\n";
} else {
	print "VALUE NOT IN DATABASE\n";
}
