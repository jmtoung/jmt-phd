#!/usr/bin/perl

use Database;
use File::Basename;

### calculates pileup on all the bam files in the database given 
### pileup is calculated only if bam file exists and a folder ending in .pileup doesn't exist

my $HOME = "/home/jmtoung/Lab";

my $ALNDB = Database->new("$HOME/aln/alnDB/alnDB_uniq.txt");
my @ALNDB = $ALNDB->read_entire_DB;

foreach my $entry (@ALNDB) {

	my $bam = $entry->[1];
	my $index = $entry->[3];

	if ($index eq '/ifs/h/toung/database/human_b36_female_gencode_refseq_ebv_50bp/human_b36_female_gencode_refseq_ebv_50bp bowtie') { 
		$index = "$HOME/database/human_b36_female_ebv.fa"; 
	} else { 
		$index = "$HOME/database/human_b36_male_ebv.fa"; 
	}

	my ($bam_name,$bam_dir,$bam_ext) = fileparse($bam,'\.bam');
	my $pileup = $bam_dir . $bam_name . ".pileup";

	next unless(-e $bam);
	next if(-e $pileup);	

	my $RUN_FILE = "sge_PileupByChrom-$bam_name.sh"; 
	open(FILE,'>$RUN_FILE') or die "[STDERR]: can't open $RUN_FILE: $!\n";

print <<"END";
#$ -cwd ### use current directory
#$ -S /usr/bin/perl ### program to execute script
#$ -M toung\@mail.med.upenn.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l mem_free=4G ### request memory

my \$ID = \$ENV{SGE_TASK_ID} - 1;

perl "$HOME/PileupByChrom.pl --bam $bam --index $index";

END

}
