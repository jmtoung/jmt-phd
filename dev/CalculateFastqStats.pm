package CalculateFastqStats;

require Exporter;

use strict;
use Getopt::Long;
use File::Basename;

our @ISA = qw(Exporter);
our @EXPORT = qw(calculate_fastq_stats $fastq_name $fastq_dir $fastq_ext);

### global variables to be exported
our ($READ_NAME_FILE, $SEQUENCE_FILE, $RETURN_TEMP_FILES);
our ($fastq_name, $fastq_dir, $fastq_ext);

sub calculate_fastq_stats {
	my $FASTQ = shift;
	### the following option allows us to just return the temp files, if we
	### wish to do some clean up (see FixFastqDuplicateRead.pl)
	$RETURN_TEMP_FILES = shift;

	-e $FASTQ or die "[STDERR]: $FASTQ doesn't exist: $!\n";
	($fastq_name, $fastq_dir, $fastq_ext) = fileparse($FASTQ,'\.fastq(\.gz)?');
	my $FASTQ_FH;
	if ($fastq_ext eq '.fastq.gz') {
		$FASTQ_FH = "gunzip -c $FASTQ |";
	} elsif ($fastq_ext eq '.fastq') {
		$FASTQ_FH = $FASTQ;
	} else {
		die "[STDERR]: $FASTQ doesn't end in .fastq or .fastq.gz\n";
	}

	open(FASTQ,$FASTQ_FH) or die "[STDERR]: cannot open '$FASTQ': $!\n";

	substr($fastq_dir,0,1) eq '/' or die "[STDERR]: '$FASTQ' is not absolute file name\n";
	### make temp files ####################################################
	do {
		my $RAND = int(rand(100));
		$READ_NAME_FILE = $fastq_dir . 'temp_READ_NAME_' . $fastq_name . "_" . $RAND . $fastq_ext;
		$SEQUENCE_FILE = $fastq_dir . 'temp_SEQ_NAME_' . $fastq_name . "_" . $RAND . $fastq_ext;
		$READ_NAME_FILE =~ s/\.gz//;
		$SEQUENCE_FILE =~ s/\.gz//;
	} until (!(-e $READ_NAME_FILE) && !(-e $SEQUENCE_FILE));
	open(READ_NAME,'>'.$READ_NAME_FILE) or die "[STDERR]: cannot open '$READ_NAME_FILE': $!\n";
	open(SEQUENCE,'>'.$SEQUENCE_FILE) or die "[STDERR]: cannot open '$SEQUENCE_FILE': $!\n";
	########################################################################

	### loop through fastq file and gather some stats ######################
	my $READ_COUNT;
	my %READ_LENGTH;
	while(my $line1 = <FASTQ>) {
		defined (my $line2 = <FASTQ>) or die "[STDERR]: line2 (sequence) undefined\n";
		defined (my $line3 = <FASTQ>) or die "[STDERR]: line3 (read_name) undefined\n";
		defined (my $line4 = <FASTQ>) or die "[STDERR]: line4 (qual scores) undefined\n";
		chomp ($line1, $line2, $line3, $line4);

		print READ_NAME $line1, "\n";
		print SEQUENCE $line2, "\n";
	
		$READ_COUNT++;
		$READ_LENGTH{length($line2)}++;
	}
	########################################################################
	
	### return temp files if $RETURN_TEMP_FILES specify so #################
	return $READ_NAME_FILE, $SEQUENCE_FILE if ($RETURN_TEMP_FILES);
	########################################################################

	### calculate $UNIQUE_SEQ_COUNT ########################################
	chomp (my $UNIQUE_SEQ_COUNT = `sort $SEQUENCE_FILE | uniq -c | wc -l`);
	########################################################################

	### calculate $UNIQUE_READ_NAME_COUNT ##################################
	chomp (my $UNIQUE_READ_NAME_COUNT = `sort $READ_NAME_FILE | uniq -c | wc -l`);
	########################################################################

	### calculate $AVERAGE_READ_LENGTH, $READ_LENGTH_DISTRIBUTION ##########
	my @READ_LENGTH_DISTRIBUTION;
	my ($AVERAGE_NUM, $AVERAGE_DEN);
	foreach my $read_length (sort {$a <=> $b} keys %READ_LENGTH) {
		push(@READ_LENGTH_DISTRIBUTION,$read_length.':'.$READ_LENGTH{$read_length});
		$AVERAGE_NUM += $read_length * $READ_LENGTH{$read_length};
		$AVERAGE_DEN += $READ_LENGTH{$read_length};
	}
	my $READ_LENGTH_DISTRIBUTION = join(",",@READ_LENGTH_DISTRIBUTION);
	my $AVERAGE_READ_LENGTH = $AVERAGE_NUM / $AVERAGE_DEN;
	########################################################################

	return 	$READ_COUNT, $UNIQUE_SEQ_COUNT, $UNIQUE_READ_NAME_COUNT, 
		$AVERAGE_READ_LENGTH, $READ_LENGTH_DISTRIBUTION;
}

END {
	### delete temporary files
	unless($RETURN_TEMP_FILES) {	
		system("rm $SEQUENCE_FILE") if ($SEQUENCE_FILE && -e $SEQUENCE_FILE);	
		system("rm $READ_NAME_FILE") if ($READ_NAME_FILE && -e $READ_NAME_FILE);
	}
}

1;
