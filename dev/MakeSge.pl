#!/usr/bin/perl

use lib '/gpfs/fs121/h/toung/dev', '/home/jmtoung/Lab/dev';
use strict;
use Database;
use File::Basename;
use Getopt::Long;
use List::Util qw[min max];
use GetLineCount;

my $HOME = "/gpfs/fs121/h/toung";

my $num_files = scalar(@FASTQ);

print STDOUT <<"END";
#\$ -cwd ### use current directory
#\$ -S /usr/bin/perl ### program to execute script
#\$ -M toung\@mail.med.upenn.edu ### email address
##\$ -m ea ### mail is to be sent at abort and end time
#\$ -j y ### combine stdout and stderr
#\$ -t 1-$num_files ### array job #
#\$ -V ### use current environment variables
##\$ -pe DJ 12 ### parallel threads
#\$ -l h_vmem=6G ### request memory

use lib "$HOME/dev";
use GetLineCount;

my \$ID = \$ENV{SGE_TASK_ID} - 1;

my \$HOME = "$HOME";


print "completed\\t\$ID\\n";

END


