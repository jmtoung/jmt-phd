package CalculateRDDStats_v20121203;

require Exporter;
use strict;
use Getopt::Long;
use List::Util qw[min max sum];

our @ISA = qw(Exporter);
our @EXPORT = qw(makeRddObject getRddStats getBases makeFastq getBaseCounts calculateNumMatches);

our ($BAM,$CHROM,$POSITION,$MINQUAL,$COMBINE_TRIM,%OBJ,$i,$EDIT_DIST_FILTER,$FULL_BAM,$ALN_ID);

sub makeRddObject {
	($BAM,$CHROM,$POSITION,$MINQUAL,$COMBINE_TRIM,$EDIT_DIST_FILTER,$FULL_BAM,$ALN_ID) = @_;

	my $REGION = $CHROM . ":" . $POSITION . "-" . $POSITION;

	%OBJ = ();

	for ($i = 0; $i < @{$BAM}; $i++) {
		$BAM->[$i]->fast_pileup($REGION,\&pileupFunction);
	}

	return \%OBJ;
}

sub pileupFunction {
	my ($SEQID,$POS,$PILEUP,$SAM) = @_;

	return unless $CHROM eq $SEQID && $POS == $POSITION;
	
	foreach my $P (@{$PILEUP}) {
		my $ALN = $P->alignment;
		my $START = $ALN->pos + 1;
		my $NAME = $ALN->qname;
		my $STRAND = convertStrandName($ALN->strand);

		my $FLAG = $ALN->flag;	
		my $MAPQUAL= $ALN->qual;
		
		my $CIGAR = $ALN->cigar_str;
		my $SEQUENCE = $ALN->qseq;
		my $QUAL = $ALN->qscore; ### an array reference
	
		if (defined $ALN->get_tag_values('SP')) {
			$START = $ALN->get_tag_values('ZR');
			$CIGAR = $ALN->get_tag_values('ZC');
			$SEQUENCE = $ALN->get_tag_values('ZS');
			$QUAL = convertToQscore($ALN->get_tag_values('ZQ'));
		} 
	
		my ($POS_IN_SEQ,$BASE_LENGTH) = calculatePositionInSequence($POSITION,$START,$CIGAR);

		my $BASE = calculateBase($SEQUENCE,$POS_IN_SEQ,$BASE_LENGTH,$STRAND);
		next if $BASE eq 'S';

		my $BASE_QUAL = $QUAL->[$POS_IN_SEQ];
		next if $BASE_QUAL < $MINQUAL;

		my $PAIR = getPair($FLAG);
		my ($ADAPTER_LENGTH,$TRIMQUAL_LENGTH) = calculateTrimLength($NAME,$PAIR,$COMBINE_TRIM);

		my ($POS5,$POS3) = calculatePositionInRead($POS_IN_SEQ,$STRAND,$TRIMQUAL_LENGTH,length($SEQUENCE));
	
		my $KEY = join("|",$START,$SEQUENCE,$CIGAR); ### this essentially groups all reads by sequence and location to remove duplicate sequences
		my $ADAPTER = 0;
		if ($ADAPTER_LENGTH > 0) { $ADAPTER = 1; }

		my $SP = my $INS = my $DEL = my $MATCH = my $CLIP = 0;

		$SP = 1 if $CIGAR =~ /N/;
		$INS = 1 if $CIGAR =~ /I/;
		$DEL = 1 if $CIGAR =~ /D/;
		$CLIP = 1 if $CIGAR =~ /S/;
		$MATCH = 1 if $CIGAR =~ /^[0-9]+M$/; 
		
		my $NUM_MATCHES = calculateNumMatches($ALN->get_tag_values('MD'));

		my $EDIT_DIST;
		if ($NUM_MATCHES == -1) { ### NUM_MATCHES == -1 when MD is not defined (like from BLAT alignment)
			$EDIT_DIST = -1;
		} else {
			$EDIT_DIST = length($SEQUENCE) - $NUM_MATCHES;
			$EDIT_DIST >= 0 or die "[STDERR]: edit distance < 0 for: ", $NAME, "\n";
		}

		if (defined $EDIT_DIST_FILTER && $EDIT_DIST_FILTER) {
			my $EDIT_DIST_ALLOWED = calculateGsnapMismatchesAllowed(length($SEQUENCE));
			defined $EDIT_DIST_ALLOWED && $EDIT_DIST_ALLOWED =~ /^[0-9]+$/ or die "[STDERR]: editDistallowed not defined for $KEY\n";
											
			next if $EDIT_DIST > $EDIT_DIST_ALLOWED;
		}
				
		### VARIABLES THAT ARE BEING COUNTED
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'count'}++;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'splice'} += $SP;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'insertion'} += $INS;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'deletion'} += $DEL;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'match'} += $MATCH;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'clip'} += $CLIP;
	
		### VARIABLES WHERE YOU EXPECT DIFFERENT VALUES
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'qual'}},$BASE_QUAL);
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'mapqual'}},$MAPQUAL);
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'pos5'}},$POS5);
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'pos3'}},$POS3);
		
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'fullQual'}},convertToAscii($QUAL));
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'name'}},$NAME);
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'bam'}},$FULL_BAM->[$i]) if defined $FULL_BAM;
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'alnID'}},$ALN_ID->[$i]) if defined $ALN_ID;
		push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'posInSeq'}},$POS_IN_SEQ);
		
		### VARIABLES WHERE YOU EXPECT ONLY ONE VALUE (HENCE TALLY UP DISTRIBUTION AND CHECK LATER THAT SCALAR(@KEYS) == 1)
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'edist'}{$EDIT_DIST}++;
		$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'numMatch'}{$NUM_MATCHES}++;
	}
}

sub getBases {
	my ($OBJ,$STRANDS,$discount_indels) = @_;
		
	my %BASES;
	foreach my $STRAND (@{$STRANDS}) {
		foreach my $BASE (sort keys %{$OBJ->{$STRAND}}) {
			if (defined $discount_indels && $discount_indels) {
				next unless $BASE =~ /^[A|C|G|T]$/;
			}
			$BASE =~ tr/ACGT/TGCA/ if ($STRAND eq '-' && scalar(@{$STRANDS}) == 2);
			$BASES{$BASE}++;
		}
	}
	my @BASES = sort keys %BASES;
	return \@BASES;
}

sub getBaseCounts {
	my ($OBJ,$STRANDS,$BASES,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY,$EDIT_DIST_FILTER) = @_;
	
	my $COUNT = 0;
	foreach my $STRAND (@{$STRANDS}) {
		foreach my $B (@{$BASES}) {
			my $BASE = $B;
			$BASE =~ tr/ACGT/TGCA/ if ($STRAND eq '-' && scalar(@{$STRANDS}) == 2);
		
			next unless keys %{$OBJ->{$STRAND}{$BASE}} > 0;
		
			foreach my $ADAPTER (0,1) {
				next if $ADAPTER_ONLY && $ADAPTER == 0;
			
				foreach my $KEY (keys %{$OBJ->{$STRAND}{$BASE}{$ADAPTER}}) {
					my $count = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'count'};
				
					$COUNT++ if $UNIQUE_SEQ_ONLY;
					$COUNT += $count if !$UNIQUE_SEQ_ONLY;
				}
			}
	
		}
	}
	return $COUNT;
}

sub calculateGsnapMismatchesAllowed {
	my $length = shift;

	my $kmer = 12;

	return max(int(($length + 2)/$kmer - 2),0);
}

sub makeFastq {
	my ($OBJ,$STRANDS,$RDD_BASE,$ADAPTER_ONLY,$RDD_ID) = @_;
	
	my $READ_ID = 1;
	foreach my $STRAND (@{$STRANDS}) {
		
		my @BASES = sort keys %{$OBJ->{$STRAND}};
		
		foreach my $B (@BASES) {
		
			my $BASE = $B;
			$BASE =~ tr/ACGT/TGCA/ if ($STRAND eq '-' && scalar(@{$STRANDS}) == 2);
			
			next unless keys %{$OBJ->{$STRAND}{$BASE}} > 0;

			my $BASE_TYPE = 'ref';
			
			$BASE_TYPE = 'rdd' if $BASE eq $RDD_BASE;
			
			foreach my $ADAPTER (0,1) {
				next if $ADAPTER_ONLY && $ADAPTER == 0;

				foreach my $KEY (keys %{$OBJ->{$STRAND}{$BASE}{$ADAPTER}}) {
				
					my ($START,$SEQUENCE,$CIGAR) = split('\|',$KEY);
					
					my $SEQUENCE_ORIGINAL = $SEQUENCE; ### just so we have original copy to look at. doesn't get printed anywhere
					$SEQUENCE =~ tr/ACGT/TGCA/ if ($STRAND eq '-');
					$SEQUENCE = reverse($SEQUENCE) if ($STRAND eq '-');
					
					my $COUNT = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'count'};

					my @NUM_MATCH = keys %{$OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'numMatch'}};
					@NUM_MATCH == 1 or die "[STDERR]: multiple num_matches at $KEY\n";
					
					my $QUAL = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'fullQual'};
					my $NAME = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'name'};
					my $alnID = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'alnID'};
					my $posInSeq = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'posInSeq'};

					scalar(@{$QUAL}) == scalar(@{$NAME}) or die "[STDERR]: diff # of entries in QUAL/NAME at $CHROM:$POSITION\n";
					scalar(@{$QUAL}) == scalar(@{$alnID}) or die "[STDERR]: diff # of entries in QUAL/alnID at $CHROM:$POSITION\n";
					scalar(@{$QUAL}) == scalar(@{$posInSeq}) or die "[STDERR]: diff # of entries in QUAL/posInSeq at $CHROM:$POSITION\n";
					
					for (my $i = 0; $i < @{$QUAL}; $i++) {
						my @ENTRY;					

						my $posInSeqConv = $posInSeq->[$i];
						$posInSeqConv = length($SEQUENCE) - $posInSeqConv - 1 if $STRAND eq '-';
						my $ENTRY_NAME = '@' . join(";",$NAME->[$i],$alnID->[$i],$RDD_ID,$READ_ID,$BASE_TYPE,$posInSeqConv,$NUM_MATCH[0]);
						push(@ENTRY,$ENTRY_NAME);
						
						push(@ENTRY,$SEQUENCE);
						push(@ENTRY,'+');

						my $ENTRY_QUAL = $QUAL->[$i];
						$ENTRY_QUAL = reverse($ENTRY_QUAL) if $STRAND eq '-';
						push(@ENTRY,$ENTRY_QUAL);

						print join("\n",@ENTRY), "\n"; 
					}
					
					$READ_ID++;
				}
			}
		}
	}
}

sub getRddStats {
	my ($OBJ,$BASES,$STRANDS,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY) = @_;

	my @GROUP1 = ('count','splice','insertion','deletion','match','clip');
	my @GROUP2 = ('qual','mapqual','pos5','pos3','edist','numMatch');
	
	my %RESULTS;
	foreach my $G1 (@GROUP1) { $RESULTS{$G1} = 0; }
	foreach my $G2 (@GROUP2) { my @BLANK; $RESULTS{$G2} = \@BLANK; }
		
	foreach my $STRAND (@{$STRANDS}) {

		next unless keys %{$OBJ->{$STRAND}} > 0;

		foreach my $B (@{$BASES}) {

			my $BASE = $B;
			$BASE =~ tr/ACGT/TGCA/ if ($STRAND eq '-' && scalar(@{$STRANDS}) == 2);
			next unless keys %{$OBJ->{$STRAND}{$BASE}} > 0;

			foreach my $ADAPTER (0,1) {

				next if $ADAPTER_ONLY && $ADAPTER == 0;

				foreach my $KEY (keys %{$OBJ->{$STRAND}{$BASE}{$ADAPTER}}) {

					foreach my $G1 (@GROUP1) {

						my $VAL = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{$G1};
						if ($UNIQUE_SEQ_ONLY) {
							$RESULTS{$G1}++ if $VAL > 0;
						} else {
							$RESULTS{$G1} += $VAL;
						}
					}

					foreach my $G2 (@GROUP2) {
						my $VAL = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{$G2};
						if (ref($VAL) eq 'HASH') { 
							my @KEYS = sort {$a<=>$b} keys %{$VAL};
							scalar(@KEYS) == 1 or die "[STDERR]: more than one key for $KEY for $G2\n" if ($G2 eq 'edist' || $G2 eq 'numMatch');
							$VAL = \@KEYS;
						}
						
						if ($UNIQUE_SEQ_ONLY) {
							push(@{$RESULTS{$G2}},mean($VAL));
						} else {
							push(@{$RESULTS{$G2}},@{$VAL});	
						}
					}
					
				}
			}
		}
	}
	my @RETURN;
	foreach my $G (@GROUP1,@GROUP2) {
		push(@RETURN,$RESULTS{$G});
	}
	return \@RETURN;
}

sub mean {
	my $ARRAY = shift;
	return sum(@{$ARRAY})/scalar(@{$ARRAY});
}

sub calculatePositionInRead {
	my ($POS_IN_SEQ,$STRAND,$TRIMQUAL_LENGTH,$SEQUENCE_LENGTH) = @_;
	
	my ($POS5,$POS3); 

	if ($STRAND eq '+') {
		$POS5 = $POS_IN_SEQ + 1;
		$POS3 = $SEQUENCE_LENGTH - $POS_IN_SEQ + $TRIMQUAL_LENGTH;
	} elsif ($STRAND eq '-') {
		$POS3 = $POS_IN_SEQ + 1 + $TRIMQUAL_LENGTH;
		$POS5 = $SEQUENCE_LENGTH - $POS_IN_SEQ;
	} 
	
	return $POS5,$POS3;
}

sub calculatePositionInSequence {
	my ($POSITION,$START,$CIGAR) = @_;

	my @CIGAR = split('(M|N|S|D|I|H)',$CIGAR);
	my $POS_IN_SEQ = $POSITION - $START;
	
	my $BASE_LENGTH = 1;

	my @END = ($START - 1);
	while(@CIGAR) {
		my $LENGTH = shift(@CIGAR);
		my $TYPE = shift(@CIGAR);
				
		next if $TYPE eq 'H';
		if ($TYPE eq 'S') { 
			$POS_IN_SEQ += $LENGTH;
			next;
		}
		push(@END,$END[-1] + $LENGTH) unless $TYPE eq 'I'; 

		if ($POSITION <= $END[-1]) {
			if ($TYPE eq 'M') {
				$BASE_LENGTH += $CIGAR[0] if ($POSITION == $END[-1] && defined $CIGAR[1] && $CIGAR[1] eq 'I');
				return $POS_IN_SEQ,$BASE_LENGTH;
			} elsif ($TYPE eq 'D') {
				$POS_IN_SEQ = $POS_IN_SEQ - ($POSITION - $END[-2]);
				return $POS_IN_SEQ,-2;
			} elsif ($TYPE eq 'N') {
				return -1,-1;
			}
		} else { 
			$POS_IN_SEQ -= $LENGTH if ($TYPE eq 'N' or $TYPE eq 'D');
			$POS_IN_SEQ += $LENGTH if $TYPE eq 'I';
		}
	}
	die "[STDERR]: error in calculatePositionInSequence: $POSITION\n";
}

sub getStrand {
	my $FLAG = shift;
	
	return '-' if ($FLAG & 16);
	return '+';
}

sub getPair {
	my ($FLAG) = @_;
	
	return 1 if ($FLAG & 64);
	return 2 if ($FLAG & 128);
	return undef;
}

sub calculateBase {
	my ($SEQUENCE,$POS_IN_SEQ,$BASE_LENGTH,$STRAND) = @_;
	
	if ($BASE_LENGTH >= 1) {
		my $BASE = substr($SEQUENCE,$POS_IN_SEQ,$BASE_LENGTH);
		$BASE =~ tr/ACGT/TGCA/ if $STRAND eq '-';
		return $BASE;
	} elsif ($BASE_LENGTH == -1) {
		return 'S';
	} elsif ($BASE_LENGTH == -2) {
		return 'X';
	}
}

sub calculateTrimLength {
	my ($NAME,$PAIR,$COMBINE_TRIM) = @_;
	
	my @NAME = split('\|',$NAME);
	my $ADAPTER_LENGTH = my $TRIMQUAL_LENGTH = 0;

	if (@NAME == 2) {
		my @TRIM = split(':',$NAME[1]);
		if (@TRIM == 2) {
			$TRIMQUAL_LENGTH = $TRIM[$PAIR - 1];
		} elsif (@TRIM == 1) {
			$TRIMQUAL_LENGTH = $TRIM[0];
		}
	} elsif (@NAME >= 3) {
		if ($COMBINE_TRIM) {
			for (my $i = 1; $i < @NAME; $i++) {
				if ($NAME[$i] =~ /^[0-9]+$/) {
					$TRIMQUAL_LENGTH += $NAME[$i];
				} else {
					die "[STDERR]: invalid trimqual length for $NAME\n";
				}
			}
		} else {
			$ADAPTER_LENGTH = $NAME[1];
			$TRIMQUAL_LENGTH = $NAME[2];
		}
	}
	
	return ($ADAPTER_LENGTH,$TRIMQUAL_LENGTH);
}

sub calculateNumMatches {
	my $MD = shift;
	
	return -1 if !defined $MD;
	
	my @MD = split('([^0-9])+',$MD);
	
	my $NUM_MATCHES = 0;
	foreach my $md (@MD) {
		$NUM_MATCHES += $md if $md =~ /^[0-9]+$/;
	}
	return $NUM_MATCHES;
}

sub convertQuality {
	my $BASE_QUAL = shift;
	
	return ord($BASE_QUAL) - 33;
}

sub convertStrandName {
	my $STRAND = shift;
	
	return '+' if ($STRAND == 1);
	return '-' if ($STRAND == -1);
	die "[STDERR]: invalid strand $STRAND\n";

}

sub convertToQscore {
	my $BASE_QUAL = shift;
	
	my @BASE_QUAL = split('',$BASE_QUAL);
	my @RETURN = map(convertQuality($_),@BASE_QUAL);
	return \@RETURN;
}

sub convertToAscii {
	my $QUAL = shift;
	
	my @RETURN = map(convertPhredToAscii($_),@{$QUAL});
	
	return join("",@RETURN);	
}

sub convertPhredToAscii {
	my $PHRED = shift;
	
	return chr($PHRED + 33);
}

1;



