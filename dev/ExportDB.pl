#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib '/ifs/h/toung/dev';
use Database;

my $database;
my $export_file;

my $options = GetOptions(	"database=s" => \$database,
				"export_file=s" => \$export_file, 
			);

my $DB = Database->new($database);
$DB->exportdb($export_file);
