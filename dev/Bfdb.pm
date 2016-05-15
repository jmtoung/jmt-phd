package Bfdb;

require Exporter;
use strict;
use DBI;

our @ISA = qw(Exporter);
our @EXPORT = ('connect2db','prepare','execute','commit','disconnect');

our $dbh;

sub connect2db {
	my $dsn = 'dbi:Pg:dbname=bfdb;host=genomicsdb.med.upenn.edu;port=5432';
	$dbh = DBI->connect($dsn,'jonathan','gms46xy',{RaiseError => 1, AutoCommit => 0}) or die "database connection not made: $DBI::errstr";
	return $dbh;
}

sub prepare {
	my ($dbh,$command) = @_;
	
	return $dbh->prepare($command);
}

sub execute {
	my ($sth,$values) = @_;
	
	if (defined $values) {
		$sth->execute(@{$values});
	} else {
		$sth->execute();
	}

}

sub commit {
	my $dbh = shift;
	
	$dbh->commit();	
}

sub disconnect {
	my $dbh = shift;
	
	$dbh->disconnect();
}

1;
