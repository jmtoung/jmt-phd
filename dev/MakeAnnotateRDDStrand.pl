#!/usr/bin/perl

use lib '/ifs/h/toung/dev', '/home/jmtoung/Lab/dev';

use Database;
use File::Basename;
use Getopt::Long;
use File::Find;

my $HOME = "/ifs/h/toung";

my $directory;
our $annotation_file;

my $options = GetOptions (
	"directory=s" => \$directory,
	"annotation_file=s" => \$annotation_file
);

find(\&process_file,$directory);

sub process_file {
	next unless /\.rdd/;

	my $RUN_FILE = "AnnotateRDDStrand-$_.sh"; 

	open(RUN_FILE,">$RUN_FILE") or die "[STDERR]: can't open $RUN_FILE: $!\n";

print RUN_FILE <<"END";
#\$ -cwd ### use current directory
#\$ -S /bin/bash ### program to execute script
#\$ -M toung\@mail.med.upenn.edu ### email address
##\$ -m ea ### mail is to be sent at abort and end time
#\$ -j y ### combine stdout and stderr
#\$ -t 1-1 ### array job #
#\$ -V ### use current environment variables
##\$ -pe DJ 12 ### parallel threads
##\$ -l mem_free=4G ### request memory

perl $HOME/dev/AnnotateRDDStrand.pl --rdd_file $File::Find::name --annotation_file $annotation_file

END
	close(RUN_RILE);

}
