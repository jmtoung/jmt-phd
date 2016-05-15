#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib '/ifs/h/toung/dev';
use Database;

my $database;
my $entry;
my $primary_key;

my $options = GetOptions(	"database=s" => \$database, 
				"entry=s" => \$entry, 
				"primary_key=s" => \$primary_key
			);

my @entry = split(',',$entry);

my $DB = Database->new($database);
$DB->add2db(\@entry,$primary_key);
