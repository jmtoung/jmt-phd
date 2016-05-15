#!/usr/bin/perl -w

use strict;
use File::Basename;

my %data;

while(<STDIN>) {
	chomp;
	
	my $file = $_;
	my ($fileName,$fileDir,$fileExt) = fileparse($file,'');
	
	my @fileName = split('\.',$fileName);
	my $jobID = $fileName[-2];
	$jobID =~ s/^o//;
	my $taskID = $fileName[-1];
		
	open(FILE,$file) or die "[STDERR]: can't run $file: $!\n";
	my $error = 0;
	my $numCompleted = 0;
	my $fastq;
	while(<FILE>) {
		chomp;
		
		if ($.==1) {
			my @split = split('\s+');
			$fastq = $split[2];
			$fastq =~ /\.fastq(\.gz)?/ or die "[STDERR]: fastq doesn't end in fastq $fastq\n";
		}

		if ($_ =~ /^Processed\s([0-9]+)\squeries\sin\s[0-9]+\.[0-9]+\sseconds/) {
			$numCompleted = $1;
		}
		if ($_ =~ /error/ || $_ =~ /STDERR/ || $_ =~ /aborted/ || $_ =~ /core\sdumped/) {
			$error = 1;
		} 
	}
	close(FILE);

	$data{$fastq}{'completed'}{$taskID}{$jobID}{$numCompleted}++ if !$error && $numCompleted;
	$data{$fastq}{'error'}{$taskID}{$jobID}++ if $error && !$numCompleted;

}

my %completed;
foreach my $fastq (sort keys %data) {
	print "-------------------------\n";
	foreach my $type ('completed','error') {
		my @taskID = sort {$a<=>$b} keys %{$data{$fastq}{$type}};
		print join("\t",$fastq,$type,scalar(@taskID)), "\n";
		if ($type eq 'completed') {
			foreach my $taskID (@taskID) {
				$completed{$taskID}++;
				foreach my $jobID (sort {$a<=>$b} keys %{$data{$fastq}{$type}{$taskID}}) {
					foreach my $numCompleted (sort {$a<=>$b} keys %{$data{$fastq}{$type}{$taskID}{$jobID}}) {
						print join("\t",'',$taskID,$numCompleted), "\n";
					}		
				}
			}
		}
	}
}

my @completed = sort {$a<=>$b} keys %completed;
print "------------------------------------------------------------\n";
print "NON COMPLETED TASK IDS\n";
for (my $i=$completed[0]; $i <= $completed[-1]; $i++) {
	print "notCompleted\t$i\n" if !defined $completed{$i};
}

### get currently running...

##my $command = "qstat -u '*' |";

##open(COMMAND,$command) or die "[STDERR]: can't fork $command: $!\n";

##while(<COMMAND>) {
##	chomp;
	
##	next unless $_ =~ /^\s[0-9]+/;
	
##	my @split = split('\s+');

##	my $user = $split[4];
##	next unless $user =~ /toung/ || $user =~ /xiwa/ || $user =~ /jdevlin/ || $user =~ /catt/;

##	my $jobName = $split[3];
##	next unless $jobName =~ /sge_RunGsn/;
	
##	my $jobID = $split[1];
##	my $taskID = $split[10];

##}
