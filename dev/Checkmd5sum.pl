#!/usr/bin/perl -w

use strict;

my $directory = "/ifs/h/toung/rdd/1000_genomes";
my $md5sum = "/ifs/h/toung/rdd/1000_genomes/pilot_data/pilot_data.alignment.index";
my $sge_directory = "/ifs/h/toung/rdd/1000_genomes/pilot_data/sge";

### OPEN MD5SUM AND STORE ALL FILES TO BE PROCESSED INTO HASH
open(MD5SUM,$md5sum) or die "[stderr]: can't open md5sum: $!\n";
my %MD5SUM;
while(<MD5SUM>) {
        chomp;
	my @split = split('\t');
        $MD5SUM{$split[0]} = $split[1] if $split[0];
        $MD5SUM{$split[4]} = $split[5] if $split[4];
        $MD5SUM{$split[6]} = $split[7] if $split[6];
}

my @SUCCESS;
my @NOT_EXIST;
my @ERROR;
my @MD5SUM;
my $SGE_COUNTER = 0;

foreach my $FILE (keys %MD5SUM) {
        my $FULL_FILE = $directory . (substr($directory,-1,1) ne '/' && '/') . $FILE;

	### CHECK IF FILE EXISTS
        if(-e $FULL_FILE) {
		### IF FILE EXISTS, CHECK IF MD5SUM OF THAT FILE EXISTS
		my $md5sum = $FULL_FILE . ".md5sum";
		if(-e $md5sum) {
			open(md5sum,$md5sum) or die "[stderr:] can't open '$md5sum': $!\n";
			my $mvalue = <md5sum>;
			my @split = split('\s+',$mvalue);
			### IF md5sum FILE EXISTS, CHECK THAT IT'S EQUAL TO ONE IN HASH
			if ($split[0] eq $MD5SUM{$FILE}) {
				push(@SUCCESS,$FILE);
			### IT md5sum FILE EXISTS & VALUE IS WRONG, DELETE FILE			
			} else {
				push(@ERROR,$FILE);
				### DELETE THE FILE

			        my $RUN_FILE = $sge_directory . (substr($sge_directory,-1,1) ne '/' && '/') . "run_delete-$SGE_COUNTER.sh";
				$SGE_COUNTER++;	
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

rm -r $FULL_FILE

END
				close(RUN_FILE);
			}	
		### IF md5sum FILE DOESN'T EXIST, RUN MD5SUM
		} else {
			push(@MD5SUM,$FILE);
		        my $RUN_FILE = $sge_directory . (substr($sge_directory,-1,1) ne '/' && '/') . "run_md5sum-$SGE_COUNTER.sh";
			$SGE_COUNTER++;	
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

md5sum $FULL_FILE > $FULL_FILE.m5dsum

END
			close(RUN_FILE);
		}
	### IF FILE DOESN'T EXIST, PUSH TO NOT EXIST
        } else {
		push(@NOT_EXIST,$FILE);
        }
}

print "number success: ", scalar(@SUCCESS), "\n";
print "number md5sum: ", scalar(@MD5SUM), "\n";
print "number error: ", scalar(@ERROR), "\n";
print "number not exist: ", scalar(@NOT_EXIST), "\n";

