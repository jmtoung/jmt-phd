package Database;

use strict;

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;

	my $database = shift;
	$database .= "/" if substr($database,-1,1) ne "/";
	-d $database or die "[STDERR]: '$database' is not a valid directory\n";

	my $self = { 'database' => $database};
	bless $self, $class;
	return $self;
}

sub add2db {
	my $self = shift;
	my $entry = shift; ### a reference to an array
	my $primary_key = shift; ### this is the column of the primary key
	my $filename = shift; ### this is the name of the file (could be value of primary key, could be undefined)

	### check all variables in $entry are defined
	foreach my $ENTRY (@{$entry}) { defined $ENTRY or die "[STDERR]: undefined values => ", join(',',@{$entry}), "\n"; }

	my $int = 0;
	if (defined $filename) {
		my $file = $self->{'database'} . $filename;
		open(FILE,">$file") or die "[STDERR]: can't open $file: $!\n";
		print FILE join("\t",@{$entry}), "\n";
		close(FILE);
	} else  {
		### check if there already exists a file with the primary key value
		my $existing_file;
		opendir(DATABASE,$self->{'database'}) or die "[STDERR]: can't open $self->{'database'}: $!\n";
		my @files = readdir(DATABASE);
		close(DATABASE);
		
		foreach my $file (@files) {
			next if ($file eq '.' or $file eq '..');
			open(FILE,$self->{'database'}.$file) or die "[STDERR]: can't open $file: $!\n";
			my $file_entry = <FILE>;
			chomp $file_entry;
			my @split = split('\t',$file_entry);
			$existing_file = $file if ($split[$primary_key] eq $entry->[$primary_key - 1]);
		}

		my $file;		
		### if file already exists, then overwrite existing file
		if (defined $existing_file) { 
			$file = $self->{'database'} . $existing_file;
			unshift(@{$entry},$existing_file);
		### otherwise, give it the next existing file slot
		} else { 
			$file = $self->{'database'} . ++$int;
			while(-e $file) { $file = $self->{'database'} . ++$int; }
			unshift(@{$entry},$int);
		}
		open(FILE,">$file") or die "[STDERR]: can't open $file: $!\n";
		print FILE join("\t",@{$entry}), "\n";
		close(FILE);
	}
}

sub lookup {
	my $self = shift;
	my $col_value = shift; ### the column that contains the value
	my $value = shift; ### the value you want to match on
	my $col_lookup = shift; ### the column we want to look up

	opendir(DATABASE,$self->{'database'}) or die "[STDERR]: can't open $self->{'database'}: $!\n";
	my @files = readdir(DATABASE);
	close(DATABASE);

	foreach my $file (@files) {
		next if ($file eq '.' or $file eq '..');
		$file = $self->{'database'} . $file;
		open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";
		chomp (my $entry = <FILE>);
		my @split = split('\t',$entry);
		if ($split[$col_value] eq $value) { return $split[$col_lookup]; }
	}
	return undef;
}

sub exportdb {
	my $self = shift;
	my $export_file = shift;

	opendir(DATABASE,$self->{'database'}) or die "[STDERR]: can't open $self->{'database'}: $!\n";
	my @files = readdir(DATABASE);
	close(DATABASE);

	my @cat_files;
	foreach my $file (@files) {
		next if ($file eq '.' or $file eq '..');
		$file = $self->{'database'} . $file;
		push(@cat_files,$file);
	}

	my $cat_files = join(" ",@cat_files);
	system("cat $cat_files > $export_file");
}

sub cleandb {
	my $self = shift;
	my $col_file = shift; ### this is the column that contains the path to the file that must exist

	opendir(DATABASE,$self->{'database'}) or die "[STDERR]: can't open $self->{'database'}: $!\n";
	my @files = readdir(DATABASE);
	close(DATABASE);

	foreach my $file (@files) {
		next if ($file eq '.' or $file eq '..');
		$file = $self->{'database'} . $file;
		open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";
		chomp (my $entry = <FILE>);
		my @split = split('\t',$entry);
		system("rm $file") unless (-e $split[$col_file]); 
	}
}

sub readdb {
	my $self = shift;
	
	opendir(DATABASE,$self->{'database'}) or die "[STDERR]: can't open $self->{'database'}: $!\n";
	my @files = readdir(DATABASE);
	close(DATABASE);
	
	my @database;
	foreach my $file (@files) {
		next if ($file eq '.' or $file eq '..');
		$file = $self->{'database'} . $file;
		open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";
		chomp (my $entry = <FILE>);
		my @split = split('\t',$entry);
		push(@database,\@split);
	}
	return \@database;
}

1;
