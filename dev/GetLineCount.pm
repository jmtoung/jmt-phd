package GetLineCount;

require Exporter;

our @ISA = qw(Exporter);

our @EXPORT = qw(getLineCount getColumnCount runCommand getPwd hasHeader formatInteger);

sub hasHeader {
	my $file = shift;
	
	open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";
	my $firstLine = <FILE>;
	if ($firstLine =~ /^#/) {
		return 1;
	} else {
		return 0;
	}
}

sub getColumnCount {
	my $file = shift;
	
	my %count;
	open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";	
	while(<FILE>) {
		chomp;
		my @split = split('\t');
		$count{scalar(@split)}++;
	}
	
	my @counts = sort {$a<=>$b} keys %count;
	scalar(@counts) == 1 or die "[STDERR]: more than one count for $file: " . join(",",@counts) . "\n";
	
	return $counts[0];
}

sub getLineCount {
	my $file = shift;

	defined $file or die "[STDERR]: file not defined\n";	
	-e $file or die "[STDERR]: file $file doesn't exist\n";
	
	my $wc = `wc -l $file`;
	chomp $wc;
	
	my @wc = split('\s+',$wc);
	my $lineCount = $wc[0];
	
	$lineCount =~ /^[0-9]+$/ or die "[STDERR]: lineCount not a number $lineCount\n";

	return $lineCount;
	
}

sub runCommand {
	my $command = shift;
	
	print STDERR $command, "\n";
	!system($command) or die "[STDERR]: can't run $command: $!\n";
	
}

sub formatInteger {
	my ($integer,$numTotalZeros) = @_;
	defined $numTotalZeros or die "[STDERR]: numTotalZeros not defined\n";	
	return (0 x ($numTotalZeros - length($integer))) . $integer;
}

sub getPwd {
	my $pwd = `pwd`;
	chomp $pwd;
	return $pwd;
}
