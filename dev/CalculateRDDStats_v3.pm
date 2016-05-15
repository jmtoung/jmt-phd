package CalculateRDDStats_v3;

require Exporter;
use strict;
use Getopt::Long;
use List::Util qw(sum);

our @ISA = qw(Exporter);
our @EXPORT = qw(make_RDD_stats_object get_RDD_stats get_bases);

our ($BAM,$CHROM,$POSITION,$MINQUAL,$COMBINE_TRIM,%OBJ,$bam);

sub make_RDD_stats_object {
	($BAM,$CHROM,$POSITION,$MINQUAL,$COMBINE_TRIM) = @_;
	my $REGION = $CHROM . ":" . $POSITION . "-" . $POSITION;

	%OBJ = ();
	foreach $bam (@{$BAM}) {
		$bam->fast_pileup($REGION,\&callback);
	}

	return \%OBJ;
}

sub callback {
	my ($SEQID,$POS,$PILEUP,$SAM) = @_;

	return unless $CHROM eq $SEQID && $POS == $POSITION;

	foreach my $P (@{$PILEUP}) {
		my $ALN = $P->alignment;
		my $START = $ALN->pos + 1; ### add one to start position

		my $STRAND = $ALN->strand;
		if ($STRAND == 1) { $STRAND = '+'; } 
		elsif ($STRAND == -1) { $STRAND = '-'; }
		else { die "[STDERR]: invalid strand $STRAND\n"; }
		
		my $SEQUENCE = $ALN->qseq;
		my $MAPQUAL= $ALN->qual;
		my $CIGAR = $ALN->cigar_str;
		
		### the 0-based position of read in sequence and the length of the sequence (>1 for insertions)
		### if we're at deletion, $POSITION_IN_READ is the leftmost non-deleted base, so we get the quality score of that base, $BASE will be 'X'
		### if we're at splice position, $POSITION_IN_READ is whatever, $LENGTH is -1 and we skip
		### insertions are tacked onto the left base

		my ($POSITION_IN_READ,$LENGTH) = calculate_position_in_read($POSITION,$START,$CIGAR);
		
		my $BASE;
		if ($LENGTH >= 1) {
			$BASE = substr($SEQUENCE,$POSITION_IN_READ,$LENGTH);
		} elsif ($LENGTH == -1) {
			next;
		} elsif ($LENGTH == -2) {
			$BASE = 'X';
		}

		### complement sequence if on negative strand
		$BASE =~ tr/ACGT/TGCA/ if $STRAND eq '-';
		
		my $READ_QUAL = $ALN->qscore->[$POSITION_IN_READ];
		next if $READ_QUAL < $MINQUAL;

		my @READ_NAME = split('\|',$ALN->qname);
		my $ADAPTER_LENGTH = 0;
		my $TRIMQUAL_LENGTH = 0;

		if (@READ_NAME == 2) {
			my @TRIM = split(':',$READ_NAME[1]);
			if (@TRIM == 2) {
				my $PAIR;
				my $FLAG = $ALN->flag;

				if ($FLAG & 64) {
					$PAIR = 1;
				} elsif ($FLAG & 128) {
					$PAIR = 2;
				} else {
					die "[STDERR]: undefined pair for $ALN->qname (flag is $FLAG)\n";
				}
				$TRIMQUAL_LENGTH = $TRIM[$PAIR - 1];
			} elsif (@TRIM == 1) {
				$TRIMQUAL_LENGTH = $TRIM[0];
			}
		} elsif (@READ_NAME >= 3) {
			if ($COMBINE_TRIM) {
				for (my $i = 1; $i < @READ_NAME; $i++) {
					if ($READ_NAME[$i] =~ /^[0-9]+$/) {
						$TRIMQUAL_LENGTH += $READ_NAME[$i];
					} else {
						die "[STDERR]: invalid trimqual length for $ALN->qname\n";
					}
				}
			} else {
				$ADAPTER_LENGTH = $READ_NAME[1];
				$TRIMQUAL_LENGTH = $READ_NAME[2];
			}
		}

		my ($POS5,$POS3) = calculate_positions($POSITION_IN_READ,$STRAND,$TRIMQUAL_LENGTH,length($SEQUENCE));
		my $KEY = $START . "|" . $SEQUENCE . "|" . $CIGAR; ### this essentially groups all reads by sequence and location to remove duplicate sequences
		my $ADAPTER = 0;
		if ($ADAPTER_LENGTH > 0) { $ADAPTER = 1; }

		### see if it's spliced
		my $SP = my $INS = my $DEL = my $MATCH = my $CLIP = 0;

		$SP = 1 if ($CIGAR =~ /N/ or defined $ALN->get_tag_values('SP'));
		$INS = 1 if $CIGAR =~ /I/;
		$DEL = 1 if $CIGAR =~ /D/;
		$CLIP = 1 if $CIGAR =~ /S/;
		unless ($SP) { $MATCH = 1 if $CIGAR =~ /^[0-9]+M$/; }

		### get the edit distance of reads
		my $EDITDIST = get_edit_dist($ALN->get_tag_values('MD'));
		
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'count'}++;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'splice'} += $SP;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'insertion'} += $INS;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'deletion'} += $DEL;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'match'} += $MATCH;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'clip'} += $CLIP;
		
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'qual'}},$READ_QUAL);
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'mapqual'}},$MAPQUAL);
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'pos5'}},$POS5);
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'pos3'}},$POS3);
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'editdist'}},$EDITDIST);

	}
}

sub get_bases {
	my ($OBJ,$STRANDS) = @_;
		
	my %BASES;
	foreach my $STRAND (@{$STRANDS}) {
		foreach my $BASE (sort keys %{$OBJ->{$STRAND}}) {
			if ($STRAND eq '-' && scalar(@{$STRANDS}) == 2) {
				$BASE =~ tr/ACGT/TGCA/;
			}
			$BASES{$BASE}++;
		}
	}
	my @BASES = sort keys %BASES;
	return \@BASES;
}

sub get_RDD_stats {
	my ($OBJ,$BASES,$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY) = @_;

	my $COUNT = my $SPLICE = my $INSERTION = my $DELETION = my $MATCH = my $CLIP = 0;
	my (@QUAL,@MAPQUAL,@POS5,@POS3,@EDITDIST);

	foreach my $STRAND (@{$STRANDS}) {
		next unless keys %{$OBJ->{$STRAND}} > 0;
		foreach my $BASE (@{$BASES}) {

			if ($STRAND eq '-' && scalar(@{$STRANDS}) == 2) {
				$BASE =~ tr/ACGT/TGCA/;
			}
			next unless keys %{$OBJ->{$STRAND}{$BASE}} > 0;

			foreach my $ADAPTER (0,1) {

				next if $ADAPTER_ONLY && $ADAPTER == 0;

				foreach my $KEY (keys %{$OBJ->{$STRAND}{$BASE}{$ADAPTER}}) {

					my $QUAL = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'qual'};
					my $MAPQUAL = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'mapqual'};
					my $POS5 = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'pos5'};
					my $POS3 = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'pos3'};
					my $EDITDIST = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'editdist'};

					my $AVERAGE_QUAL = mean($QUAL);					
					my $AVERAGE_MAPQUAL = mean($MAPQUAL);
					my $AVERAGE_POS5 = mean($POS5);
					my $AVERAGE_POS3 = mean($POS3);
					my $AVERAGE_EDITDIST = mean($EDITDIST);
					
					if ($UNIQUE_SEQ_ONLY) {
						push(@QUAL,$AVERAGE_QUAL);
						push(@MAPQUAL,$AVERAGE_MAPQUAL);
						push(@POS5,$AVERAGE_POS5);
						push(@POS3,$AVERAGE_POS3);
						push(@EDITDIST,$AVERAGE_EDITDIST);
						$COUNT++;
						$SPLICE++ if $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'splice'} > 0;
						$INSERTION++ if $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'insertion'} > 0;
						$DELETION++ if $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'deletion'} > 0;
						$MATCH++ if $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'match'} > 0;
						$CLIP++ if $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'clip'} > 0;
					} else {
						push(@QUAL,@{$QUAL});
						push(@MAPQUAL,@{$MAPQUAL});
						push(@POS5,@{$POS5});
						push(@POS3,@{$POS3});
						push(@EDITDIST,@{$EDITDIST});
						$COUNT += $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'count'};
						$SPLICE += $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'splice'};
						$INSERTION += $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'insertion'};
						$DELETION += $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'deletion'};
						$MATCH += $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'match'};
						$CLIP += $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'clip'};
					}
				}
			}
		}
	}
	return ($COUNT,\@QUAL,\@MAPQUAL,\@POS5,\@POS3,\@EDITDIST,$SPLICE,$INSERTION,$DELETION,$MATCH,$CLIP);
}

sub get_edit_dist {
	my $MD = shift;
	
	my $EDITDIST = 0;
	my @MD = split('([^0-9])',$MD);
	foreach my $md (@MD) {
		$EDITDIST += length($md) if $md =~ /^[A|C|G|T]$/;
	}
	return $EDITDIST;
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
sub calculate_position_in_read {
	my ($POSITION,$START,$CIGAR) = @_;

	my @CIGAR = split('(M|N|S|D|I|H)',$CIGAR);
	my $POSITION_IN_READ = $POSITION - $START;
	
	my $NUM_BASES = 1; ### this is the number of bases or length of the sequence (>1 for insertions)

	my @END = ($START - 1); ### keep track of the ends of each piece of the read;
	while(@CIGAR) {
		my $LENGTH = shift(@CIGAR); ### This is the length of this piece of the read
		my $TYPE = shift(@CIGAR); ### This is the type of this piece of the read (M for exon, N for intron)
				
		next if $TYPE eq 'H'; ### added on 4/13/2012 to deal with hard clippings in DNA-Seq
		if ($TYPE eq 'S') { ### 
			$POSITION_IN_READ += $LENGTH;
			next;
		}
		push(@END,$END[-1] + $LENGTH) unless $TYPE eq 'I'; ### if type is 'I' there was no "piece"

		if ($POSITION <= $END[-1]) { ### Check if your position falls within this piece of the read
			if ($TYPE eq 'M') {
				### If the type we're on is 'M', that means read covers exon
				### check if there is a next type ($CIGAR[1]) and length ($CIGAR[0]), AND if next type is an insertion...if so, change $NUM_BASES to $CIGAR[0]
				$NUM_BASES += $CIGAR[0] if ($POSITION == $END[-1] && defined $CIGAR[0] && defined $CIGAR[1] && $CIGAR[1] eq 'I');
				return $POSITION_IN_READ,$NUM_BASES;
			} elsif ($TYPE eq 'D') {
				### if it's a deletion, return the left-most position that aligns to reference (not current end, but previous end -> $END[-2])
				$POSITION_IN_READ -= $POSITION - $END[-2];
				return $POSITION_IN_READ,-2;
			} elsif ($TYPE eq 'N') {
				### if it's a splice, return -1 and -1...
				return -1,-1;
			}
		} else { 
			### If our position is not within this piece of the read...
			### check if type is 'N', which means we're hovering over N so we need to subtract
			### check if type is 'D', which is the same thing...need to subtract gap
			$POSITION_IN_READ -= $LENGTH if ($TYPE eq 'N' or $TYPE eq 'D');
			### check if type is 'I'. if so we need to add to position in read
			$POSITION_IN_READ += $LENGTH if $TYPE eq 'I';
		}
	}
	die "[STDERR]: error calculating position in read: $POSITION\n";
}

1;
