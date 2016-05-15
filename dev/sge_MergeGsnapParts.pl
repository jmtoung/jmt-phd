#$ -cwd ### use current directory
#$ -S /usr/bin/perl ### program to execute script
#$ -M toung@mail.med.upenn.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
##$ -l mem_free=10G ### request memory

use strict;

my $ID = $ENV{SGE_TASK_ID} - 1;

my @folders = qw(
/gpfs/fs121/h/toung/aln/baseline/GM12750_onelane_trimlowqual
);

my $tag = "defv20120410";

opendir(DIR,$folders[$ID]) or die "[STDERR]: can't open directory: $folders[$ID]: $!\n";

my $output;
my %bam;
while(my $file = readdir(DIR)) {
	next unless $file =~ /^gsnap.*$tag\.[0-9]+\.bam$/;

	my @split = split('\.',$file);
	$split[2] eq $tag or die "[STDERR]: doesn't equal $tag\n";
	my $part = int($split[3]);
	
	my @output = @split;
	$output[3] = $output[4];
	pop(@output);
	$output = $folders[$ID] . "/" . join(".",@output);
	
	my $new_file = $folders[$ID] . "/" . join(".",@split);
	
	my $full_file = $folders[$ID] . "/" . $file;

	### rename from restart2 to def
	my $command = "mv $full_file $new_file";

	print "moving\t$full_file\n";
	print "to\t$new_file\n\n";
	
	!system($command) or die "[STDERR]: can't run $command: $!\n";
	
	$bam{$part} = $new_file;
}

my @bam_files;
foreach my $part (sort {$a<=>$b} keys %bam) {
	push(@bam_files,$bam{$part});
}

print "\n\nFILES_ARE...\n";
print join("\n",@bam_files), "\n";

my $bam_files = join(",",@bam_files);

my $home = "/gpfs/fs121/h/toung/oldhome";

my $command = "perl $home/dev/MergeBamFiles.pl --bam_files $bam_files | samtools view -bS - > $output";

print "\n\nCOMMAND_IS...\n";

print $command, "\n";

!system($command) or die "[STDERR]: can't run $command: $!\n";

print "completed\t$ID\n";
