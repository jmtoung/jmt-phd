#$ -cwd ### use current directory
#$ -S /usr/bin/perl ### program to execute script
#$ -M toung@mail.med.upenn.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l h_vmem=2G ### request memory

use strict;
use File::Basename;
use lib '/home/jmtoung/Lab/dev', '/gpfs/fs121/h/toung/dev';
use GetLineCount;

my $ID = $ENV{SGE_TASK_ID} - 1;
