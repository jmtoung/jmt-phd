#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib '/ifs/h/toung/dev','/home/jmtoung/Lab/dev';
use CalculateFastqStats;
use Database;
use Timestamp;

################################################################################
### This script takes in a fastq file (--fastq) and uploads the data into a data-
### base, which is contained in the fastq/fastqDB folder. The name of the fastq
### file (--fastq) is the primary key in this database; an entry will be overwri-
### tten it has the same fastq file name as the one you are inserting. The parent
### fastq file (--original_fastq), for example, in the case of a trimmed fastq file,
### the parent file is the untrimmed fastq file, must also be supplied along with
### the sample name (--sample). Various parameters, (--read_count, --unique_seq_count,...)
### can be provided. If any of them are not provided, the calculate_fastq_stats
### function in the CalculateFastqStats library will be called.
################################################################################

my $database;
my $fastq;
my $original_fastq;
my $sample;
my $read_count;
my $unique_seq_count;
my $read_name_count;
my $unique_read_name_count;
my $average_read_length;
my $read_length_distribution;

my $options = GetOptions(
	"database=s" => \$database,
	"fastq=s" => \$fastq, 
	"original_fastq=s" => \$original_fastq, 
	"sample=s" => \$sample, 
	"read_count=i" => \$read_count, 
	"unique_seq_count=i" => \$unique_seq_count,			
	"unique_read_name_count=i" => \$unique_read_name_count,
	"average_read_length=i" => \$average_read_length,
	"read_length_distribution=s" => \$read_length_distribution,
);

### IF ANY OF FOLLOWING !DEFINED, CALL FUNCTION TO DEFINE ALL OF THEM ##########
if (	!defined $read_count || !defined $unique_seq_count ||
	!defined $unique_read_name_count || !defined $average_read_length ||
	!defined $read_length_distribution) {
	($read_count,$unique_seq_count,$unique_read_name_count,
	$average_read_length,$read_length_distribution) = calculate_fastq_stats($fastq);
}
################################################################################

################################################################################
my $timestamp = get_timestamp;
my @entry = ($fastq,$original_fastq,$sample,$read_count,$unique_seq_count,
	$unique_read_name_count,$average_read_length,$read_length_distribution,
	$timestamp);
################################################################################

### DATABASE PARAMETERS -- UPDATE THIS #########################################
my $fastqDB = Database->new($database);
$fastqDB->add_entry(\@entry,1);
################################################################################
