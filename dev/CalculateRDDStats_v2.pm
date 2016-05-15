package CalculateRDDStats_v2;

require Exporter;
use lib '/home/jmtoung/Lab/dev', '/ifs/h/toung/dev', '/gpfs/fs121/h/toung/oldhome/dev';
use strict;
use Getopt::Long;
use ComplementBase;
use List::Util qw(sum);

our @ISA = qw(Exporter);
our @EXPORT = qw(make_RDD_stats_object get_RDD_stats get_bases);

our ($BAM,$CHROM,$POSITION,$MINQUAL,%OBJ);

### this method generates an "object"
sub make_RDD_stats_object {
	($BAM,$CHROM,$POSITION,$MINQUAL) = @_;
	my $REGION = $CHROM . ":" . $POSITION . "-" . $POSITION;

	%OBJ = ();
	$BAM->fast_pileup($REGION,\&callback);

	return \%OBJ;
}

sub callback {
	my ($SEQID,$POS,$PILEUP,$SAM) = @_;

	return unless $CHROM eq $SEQID && $POS == $POSITION; ### return unless we're on the position we want!

	foreach my $P (@{$PILEUP}) {
		my $ALN = $P->alignment;
		my $START = $ALN->pos + 1;
		my $STRAND = $ALN->strand;
		if ($STRAND == 1) { $STRAND = '+'; } 
		elsif ($STRAND == -1) { $STRAND = '-'; }
		else { die "[STDERR]: invalid strand $STRAND\n"; }
		
		my $SEQUENCE = $ALN->qseq;
		my $MAPQUAL= $ALN->qual;
		my $CIGAR = $ALN->cigar_str;

		### the 0-based position of read in sequence
		my $POSITION_IN_READ = calculate_position_in_read($POSITION,$START,$CIGAR); 
		### -1 means read base is empty
		next if $POSITION_IN_READ == -1; 
		
		my $BASE = substr($SEQUENCE,$POSITION_IN_READ,1);
		$BASE = complement_base($BASE) if $STRAND eq '-';
		my $READ_QUAL = $ALN->qscore->[$POSITION_IN_READ];
		next if $READ_QUAL < $MINQUAL;

		my @READ_NAME = split('\|',$ALN->qname);
		my $ADAPTER_LENGTH = 0;
		my $TRIMQUAL_LENGTH = 0;
		if (scalar(@READ_NAME) == 3) {
			$ADAPTER_LENGTH = $READ_NAME[1];
			$TRIMQUAL_LENGTH = $READ_NAME[2];
		} elsif (scalar(@READ_NAME) == 2) {
			$TRIMQUAL_LENGTH = $READ_NAME[1];
		}

		my ($POS5,$POS3) = calculate_positions($POSITION_IN_READ,$STRAND,$TRIMQUAL_LENGTH,length($SEQUENCE));
		my $KEY = $START . "|" . $SEQUENCE . "|" . $CIGAR; ### this essentially groups all reads by sequence and location to remove duplicate sequences
		my $ADAPTER = 0;
		if ($ADAPTER_LENGTH > 0) { $ADAPTER = 1; }

		### see if it's spliced
		my $SP = my $INS = my $DEL = my $MATCH = 0;

		$SP = 1 if defined $ALN->get_tag_values('SP');		
		$INS = 1 if $CIGAR =~ /I/;
		$DEL = 1 if $CIGAR =~ /D/;
		$MATCH = 0;
		unless ($SP) {
			$MATCH = 1 if $CIGAR =~ /^[0-9]+M$/;
		}
		
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'count'}++;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'splice'} += $SP;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'insertion'} += $INS;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'deletion'} += $DEL;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'match'} += $MATCH;
		
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'pos3'}},$POS3);
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'pos5'}},$POS5);
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'mapqual'}},$MAPQUAL);
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'qual'}},$READ_QUAL);
	}
}

sub get_bases {
	my ($OBJ,$STRANDS) = @_;
	
	my %BASES;
	foreach my $STRAND (@{$STRANDS}) {
		foreach my $BASE (keys %{$OBJ->{$STRAND}}) {
			$BASE = complement_base($BASE) if ($STRAND eq '-' && scalar(@{$STRANDS}) == 2);
			$BASES{$BASE}++;
		}
	}
	my @BASES = keys %BASES;
	return \@BASES;
}
sub get_RDD_stats {
	my ($OBJ,$BASES,$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY) = @_;

	my $COUNT = my $SPLICE = my $INSERTION = my $DELETION = my $MATCH = 0;
	my (@QUAL,@MAPQUAL,@POS5,@POS3);

	foreach my $STRAND (@{$STRANDS}) {
		next unless exists $OBJ->{$STRAND};
		foreach my $BASE (@{$BASES}) {
			
			$BASE = complement_base($BASE) if ($STRAND eq '-' && scalar(@{$STRANDS}) == 2);
			next unless exists $OBJ->{$STRAND}{$BASE};

			foreach my $ADAPTER (0,1) {

				next if $ADAPTER_ONLY && $ADAPTER == 0;

				foreach my $KEY (keys %{$OBJ->{$STRAND}{$BASE}{$ADAPTER}}) {

					my $MAPQUAL = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'mapqual'};
					my $POS5 = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'pos5'};
					my $POS3 = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'pos3'};
					my $QUAL = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'qual'};

					my $AVERAGE_MAPQUAL = mean($MAPQUAL);
					my $AVERAGE_POS5 = mean($POS5);
					my $AVERAGE_POS3 = mean($POS3);
					my $AVERAGE_QUAL = mean($QUAL);
					
					if ($UNIQUE_SEQ_ONLY) {
						push(@POS5,$AVERAGE_POS5);
						push(@POS3,$AVERAGE_POS3);
						push(@MAPQUAL,$AVERAGE_MAPQUAL);
						push(@QUAL,$AVERAGE_QUAL);
						$COUNT++;
						$SPLICE++ if $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'splice'} > 0;
						$INSERTION++ if $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'insertion'} > 0;
						$DELETION++ if $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'deletion'} > 0;
						$MATCH++ if $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'match'} > 0;
					} else {
						push(@POS5,@{$POS5});
						push(@POS3,@{$POS3});
						push(@MAPQUAL,@{$MAPQUAL});
						push(@QUAL,@{$QUAL});
						$COUNT += $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'count'};
						$SPLICE += $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'splice'};
						$INSERTION += $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'insertion'};
						$DELETION += $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'deletion'};
						$MATCH += $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'match'};
					}
				}
			}
		}
	}
	return ($COUNT,\@QUAL,\@MAPQUAL,\@POS5,\@POS3,$SPLICE,$INSERTION,$DELETION,$MATCH);
}


sub mean {
	my $ARRAY = shift;
	return sum(@{$ARRAY})/scalar(@{$ARRAY});
}

### This method takes in the position in the read calculated from the method 'calculate_read_base' along with the strand 
### to determine positions relative to 3' and 5' ends
sub calculate_positions {
	my ($POSITION_IN_READ,$STRAND,$TRIMQUAL_LENGTH,$SEQUENCE_LENGTH) = @_;
	### position of read from 5' and 3' end
	my ($POS5,$POS3); 

	if ($STRAND eq '+') {
		$POS5 = $POSITION_IN_READ + 1;
		$POS3 = $SEQUENCE_LENGTH - $POSITION_IN_READ + $TRIMQUAL_LENGTH;
	} elsif ($STRAND eq '-') {
		$POS3 = $POSITION_IN_READ + 1 + $TRIMQUAL_LENGTH;
		$POS5 = $SEQUENCE_LENGTH - $POSITION_IN_READ;
	} else { die "[STDERR]: invalid strand $STRAND\n"; }
	return $POS5,$POS3;
}

### This method returns the actual base in the read along with the position of the base in the read (0-based)
### Extra caution must be taken to deal with spliced reads
sub calculate_position_in_read {
	my ($POSITION,$START,$CIGAR) = @_;

	my @CIGAR = split('(M|N|S|D|I)',$CIGAR);
	my $POSITION_IN_READ = $POSITION - $START;

	my $END = $START - 1; ### Use this variable to calculate the end of this piece of the read. Initialize to -1.
	while(@CIGAR) {
		my $LENGTH = shift(@CIGAR); ### This is the length of this piece of the read
		my $TYPE = shift(@CIGAR); ### This is the type of this piece of the read (M for exon, N for intron)
		if ($TYPE eq 'S') {
			$POSITION_IN_READ += $LENGTH;
			next;
		}
		$END += $LENGTH unless $TYPE eq 'I'; ### Update the end of this piece of the read

		if ($POSITION <= $END) { ### Check if your position falls within this piece of the read
			### If the type we're on is 'M', that means read covers exon
			return $POSITION_IN_READ if $TYPE eq 'M';
			### If not, we're hovering over intron or deletion, so return -1 for empty
			return -1;
		} else { ### If our position is not within this piece of the read, subtract from the position in the read the length of this entire piece of the read if it's a 'N'.			
			$POSITION_IN_READ -= $LENGTH if ($TYPE eq 'N' or $TYPE eq 'D');
			$POSITION_IN_READ += $LENGTH if $TYPE eq 'I';
		}
	}
	die "[STDERR]: error calculating position in read: $POSITION\n";
}

1;
