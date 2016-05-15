#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use lib '/ifs/h/toung/dev';
use CalculateFastqStats;

my $fastq;
my $tag;

my $options = GetOptions(
	"fastq=s" => \$fastq,
	"tag=s" => \$tag
);

my ($read_name_file, $sequence_file) = calculate_fastq_stats($fastq,1);

### sort and unique the "$read_name_file" to get a list of duplicate read names
my %DUPLICATE_READS;
open(READ_NAME,"sort $read_name_file | uniq -c |");
while(<READ_NAME>) {
	chomp;
	my @split = split('\s+');
	$DUPLICATE_READS{$split[2]}++ if $split[1] != 1;
}
close(READ_NAME);
my $num_dup_reads = scalar(keys %DUPLICATE_READS);
print "number of duplicated reads for $fastq is $num_dup_reads\n";

### open the fastq file again (don't need to do any checks b/c was done in the 
### 'calculate_fastq_stats' method already) and output new file ################
my $fastq_new = $fastq_dir . $fastq_name . "_" . $tag . $fastq_ext;
$fastq_new =~ s/\.gz//;
my $FASTQ_FH;
if ($fastq_ext eq '.fastq.gz') {
	$FASTQ_FH = "gunzip -c $FASTQ |";
} elsif ($fastq_ext eq '.fastq') {
	$FASTQ_FH = $FASTQ;
} else {
	die "[STDERR]: $FASTQ doesn't end in .fastq or .fastq.gz\n";
}
open(FASTQ_NEW,'>'.$fastq_new);
select(FASTQ_NEW);
open(FASTQ,$FASTQ_FH);
### store sequences of duplicated reads.
my %DUPLICATE_READS_SEQUENCES; 
while(my $line1 = <FASTQ>) {
	my $line2 = <FASTQ>;
	my $line3 = <FASTQ>;
	my $line4 = <FASTQ>;
	chomp ($line1, $line2, $line3, $line4);

	### first check if the read is a duplicated read. if so, then...
	if (exists $DUPLICATE_READS{$line1}) {
		### next check if the sequence is different
		### two reads might have same name but different sequences
		### if reads have same sequence, then skip
		next if (exists $DUPLICATE_READS_SEQUENCES{$line1}{$line2});
		### if it doesn't, then come up with a new read name for it
		my $num = 1;
		my $line1_new;
		my @line1 = split('\/',$line1);
		do {
			my @line1_new = @line1;
			$line1_new[0] .= "-" . $num++;
			$line1_new = join("/",@line1_new);
		} until (!exists $DUPLICATE_READS{$line1_new});
		### put the new read_name in the duplicated reads hash
		$DUPLICATE_READS{$line1_new}++;
		### put the sequence under the old read name
		$DUPLICATE_READS_SEQUENCES{$line1}{$line2}++;
		### print read
		print $line1_new, "\n", $line2, "\n", $line3, "\n", $line4, "\n";
	### if read is a not a duplicated read, then just print
	} else {
		print $line1, "\n", $line2, "\n", $line3, "\n", $line4, "\n";
	}
}

END {
	### delete temporary files
	system("rm $read_name_file");
	system("rm $sequence_file");
}
