package CalculateRDDStats_v4;

require Exporter;
use strict;
use Getopt::Long;
use List::Util qw(sum);

our @ISA = qw(Exporter);
our @EXPORT = qw(makeRddObject getRddStats getBases makeFasta);

our ($BAM,$CHROM,$POSITION,$MINQUAL,$COMBINE_TRIM,%OBJ);

sub makeRddObject {
	($BAM,$CHROM,$POSITION,$MINQUAL,$COMBINE_TRIM) = @_;
	my $REGION = $CHROM . ":" . $POSITION . "-" . $POSITION;

	%OBJ = ();
	foreach my $bam (@{$BAM}) {
		my $command = "samtools view $bam $REGION |";
		open(COMMAND,$command) or die "[STDERR]: can't open $command: $!\n";
		while(<COMMAND>) {
			chomp;
			my @split = split('\t');
			callback(\@split);
		}
	}

	return \%OBJ;
}

sub callback {
	my $split = shift;

	return unless $CHROM eq $split->[2];
	
	my $NAME = $split->[0];
	my $FLAG = $split->[1];
	my $START = $split->[3];

	my $STRAND = getStrand($FLAG);
		
	my $MAPQUAL= $split->[4];
	my $CIGAR = $split->[5];
	my $SEQUENCE = $split->[9];
	my $QUAL = $split->[10];
	
	my ($POS_IN_SEQ,$BASE_LENGTH) = calculatePositionInSequence($POSITION,$START,$CIGAR);
		
	my $BASE = calculateBase($SEQUENCE,$POS_IN_SEQ,$BASE_LENGTH,$STRAND);
	return if $BASE eq 'S';

	my $BASE_QUAL = convertQuality(substr($QUAL,$POS_IN_SEQ,1));
	return if $BASE_QUAL < $MINQUAL;

	my $PAIR = getPair($FLAG,$split);
	my ($ADAPTER_LENGTH,$TRIMQUAL_LENGTH) = getTrimLength($NAME,$PAIR,$COMBINE_TRIM);

	my ($POS5,$POS3) = calculatePositionInRead($POS_IN_SEQ,$STRAND,$TRIMQUAL_LENGTH,length($SEQUENCE));
	
	my $KEY = join("|",$START,$SEQUENCE,$STRAND,$CIGAR); ### this essentially groups all reads by sequence and location to remove duplicate sequences
	my $ADAPTER = 0;
	if ($ADAPTER_LENGTH > 0) { $ADAPTER = 1; }

	my $SP = my $INS = my $DEL = my $MATCH = my $CLIP = 0;

	$SP = 1 if $CIGAR =~ /N/;
	$INS = 1 if $CIGAR =~ /I/;
	$DEL = 1 if $CIGAR =~ /D/;
	$CLIP = 1 if $CIGAR =~ /S/;
	$MATCH = 1 if $CIGAR =~ /^[0-9]+M$/; 

	my $NUM_MATCHES = calculateNumMatches(getTagValue($split,'MD'));
	
	my $EDIT_DIST = length($SEQUENCE) - $NUM_MATCHES;

	$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'count'}++;
	$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'splice'} += $SP;
	$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'insertion'} += $INS;
	$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'deletion'} += $DEL;
	$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'match'} += $MATCH;
	$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'clip'} += $CLIP;
		
	push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'qual'}},$BASE_QUAL);
	push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'mapqual'}},$MAPQUAL);
	push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'pos5'}},$POS5);
	push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'pos3'}},$POS3);
	push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'edist'}},$EDIT_DIST);
	push(@{$OBJ{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'numMatch'}},$NUM_MATCHES);
}

sub getBases {
	my ($OBJ,$STRANDS) = @_;
		
	my %BASES;
	foreach my $STRAND (@{$STRANDS}) {
		foreach my $BASE (sort keys %{$OBJ->{$STRAND}}) {
			$BASE =~ tr/ACGT/TGCA/ if ($STRAND eq '-' && scalar(@{$STRANDS}) == 2);
			$BASES{$BASE}++;
		}
	}
	my @BASES = sort keys %BASES;
	return \@BASES;
}

sub makeFasta {
	my ($OBJ,$FASTA,$STRANDS,$REF_BASE,$RDD_BASE,$ADAPTER_ONLY,$UNIQUE_SEQ_ONLY,$RDD_ID) = @_;
	
	my $READ_ID = 1;
	foreach my $STRAND (@{$STRANDS}) {
		
		my @BASES = sort keys %{$OBJ->{$STRAND}};
		
		my $BASE_TYPE = 'tot';
		foreach my $BASE (@BASES) {
			$BASE =~ tr/ACGT/TGCA/ if ($STRAND eq '-' && scalar(@{$STRANDS}) == 2);
			
			$BASE_TYPE = 'rdd' if $BASE eq $RDD_BASE;
			
			foreach my $ADAPTER (0,1) {
				next if $ADAPTER_ONLY && $ADAPTER == 0;
				
				foreach my $KEY (keys %{$OBJ->{$STRAND}{$BASE}{$ADAPTER}}) {
				
					my ($START,$SEQUENCE,$ALN_STRAND,$CIGAR) = split('\|',$KEY);
					
					$SEQUENCE =~ tr/ACGT/TGCA/ if ($ALN_STRAND eq '-');
					$SEQUENCE = reverse($SEQUENCE);
					
					my $COUNT = $OBJ->{$STRAND}{$BASE}{$ADAPTER}{$KEY}{'count'};
					
					my $DISPLAY_ID = join(";",$COUNT,$REF_BASE,$RDD_BASE,$BASE_TYPE,$RDD_ID,$READ_ID);
					
					my $SEQ_OBJ = Bio::Seq->new(-seq => $SEQUENCE, -display_id => $DISPLAY_ID, -desc => "", -alphabet => "dna");
					
					$FASTA->write_seq($SEQ_OBJ);

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

		foreach my $BASE (@{$BASES}) {

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
	my ($FLAG,$SPLIT) = @_;
	
	return 1 if ($FLAG & 64);
	return 2 if ($FLAG & 128);
	return undef;
}

sub getTagValue {
	my ($SPLIT,$TAG) = @_;
	
	for (my $i = 11; $i < @{$SPLIT}; $i++) {
		my @TAG = split(':',$SPLIT->[$i]);
		return $TAG[2] if ($TAG[0] eq $TAG);
	}
	die "[STDERR]: tag $TAG not found in ", join("\t",@{$SPLIT}), "\n";
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

sub getTrimLength {
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
	
	my @MD = split('([^0-9])+',$MD);
	
	my $NUM_MATCHES = 0;
	foreach my $md (@MD) {
		$NUM_MATCHES += $md if $md =~ /[0-9]+/;
	}
	return $NUM_MATCHES;
}

sub convertQuality {
	my $BASE_QUAL = shift;
	
	return ord($BASE_QUAL) - 33;
}

1;



