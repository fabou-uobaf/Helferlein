#!/usr/bin/perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Data::Dumper;

## Variables
my %data = (); # holding all data

## Parameter
my $clipping_th = 3; # minimal clipping length to consider
my $initialCoverage = 1; # number of reads at the begining to call transcript attachment
my $consensusCoverage = 2; # number of reads to evaluate consensus sequence quality

GetOptions(
  "l"           => \$clipping_th,
  "c"           => \$initialCoverage,
  "cc"          => \$consensusCoverage,
  "man"         => sub{pod2usage(-verbose => 2)},
  "help|?"      => sub{pod2usage(-verbose => 1)}
    );

# read in demultiplexed, uniqued, local mapped, single-end sam/bam file
while(<>){
    next if (m/^@/);
    my @F = split/\t/, $_;

    my ($readName, $flag, $chr, $leftPos, $cigar, $readSeq)  = ($F[0], $F[1], $F[2], $F[3], $F[5], $F[9]);

    # get strand from flag
    my $strand = &parse_strand_bitflag($flag);

    # get right end from cigar
    my $rightPos = $leftPos + &get_length_from_cigar($cigar);

    # get left softclip length
    # if so: get left softclipped bases
    my ($leftClipLength, $leftClipSeq) = &get_softclip_info($cigar, $readSeq, 'L');

    # get right softclip range
    # if so: get left softclipped bases
    # if so: get left softclipped bases
    my ($rightClipLength, $rightClipSeq) = &get_softclip_info($cigar, $readSeq, 'R');

    # store in %data
    if ($leftClipLength > $clipping_th){
	$data{$chr}->{$strand}->{$leftPos}->{clippedReads}++;
	my @leftBases = split//,$leftClipSeq;
	foreach my $relativePos (0 .. $#leftBases){
	    my $effectivePos = $leftPos - $leftClipLength + $relativePos;
	    $data{$chr}->{$strand}->{$leftPos}->{clippedBases}->{$effectivePos}->{$leftBases[$relativePos]}++;
	}
    }

    if ($rightClipLength > $clipping_th){
	$data{$chr}->{$strand}->{$rightPos}->{clippedReads}++;
	my @rightBases = split//,$rightClipSeq;
	foreach my $relativePos (0 .. $#rightBases){
	    my $effectivePos = $rightPos + $relativePos;
	    $data{$chr}->{$strand}->{$rightPos}->{clippedBases}->{$effectivePos}->{$rightBases[$relativePos]}++;
	}
    }
}

 #print Dumper(%data);
# report for all polyA sites with more than $cover_th supportive reads
# the consensus sequence
# polyA length
# number of conflicting position
# number of A's
# number of positions with more than n coverage
# percentage of readbases which conflict with consensus

# print column header
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 'chr', 'position', 'strand', 'length', 'consensusSeq', 'abs_conflictPosition', 'rel_conflictPosition', 'rel_BaseCountA', 'rel_BaseCountConflicting';
foreach my $chr ( sort keys %data ){
  foreach my $strand ( sort keys %{$data{$chr}} ){
    foreach my $position ( sort {$a <=> $b} keys %{$data{$chr}->{$strand}} ){

      next if ($data{$chr}->{$strand}->{$position}->{clippedReads} < $initialCoverage); # ignore sites with less than $initialCoverage reads support

      my $length           = 0;  # attachment length
      my $conflictPos      = 0;  # positions with none perfect consensus sequence
      my $countA           = 0;  # number of seen As
      my $countTotal       = 0;  # total number of seens bases
      my $countConflicting = 0;  # number of bases in conflict with the consensus sequence (if coverage > $consensusCoverage)
      my $countNcovered    = 0;  # total number of bases (if coverage > $consensusCoverage)
      my @consensusSeq     = (); # array holding the consensus sequence

      foreach my $clipPosition ( sort {$a <=> $b} keys %{$data{$chr}->{$strand}->{$position}->{clippedBases}} ) {

        $length++;
        my $consensusBase = &hash_key_to_largest_value(\%{$data{$chr}->{$strand}->{$position}->{clippedBases}->{$clipPosition}});
        push @consensusSeq, $consensusBase;

        $conflictPos++ if (scalar (keys %{$data{$chr}->{$strand}->{$position}->{clippedBases}->{$clipPosition}}) > 1);
        $countA += $data{$chr}->{$strand}->{$position}->{clippedBases}->{$clipPosition}->{A} if (defined($data{$chr}->{$strand}->{$position}->{clippedBases}->{$clipPosition}->{A}));
        $countTotal += &hash_value_sum(\%{$data{$chr}->{$strand}->{$position}->{clippedBases}->{$clipPosition}});

        if ( $countTotal >= $consensusCoverage ) {
          foreach my $seenBase ( keys %{$data{$chr}->{$strand}->{$position}->{clippedBases}->{$clipPosition}} ) {
            $countNcovered+=$data{$chr}->{$strand}->{$position}->{clippedBases}->{$clipPosition}->{$seenBase};
            if($seenBase ne $consensusBase){
              $countConflicting+=$data{$chr}->{$strand}->{$position}->{clippedBases}->{$clipPosition}->{$seenBase};
            }
          }
        }
      }
    printf "%s\t%d\t%s\t%d\t%s\t%d\t%.2g\t%.2g\t%.2g\n", $chr, $position, $strand, $length, join("", @consensusSeq), $conflictPos, $conflictPos/$length, $countA/$countTotal, $countConflicting/$countNcovered;
    }
  }
}
###################
##  SUBROUTINES  ##
###################

sub hash_value_sum (\%) {
  my $hash = shift;
  my $value = 0;

  for my $v (values %$hash ) {
    $value += $v;
  }
  return($value);
}



sub hash_key_to_largest_value (\%) {
    my $hash = shift;
    keys %$hash;       # reset the each iterator

    my ($large_key, $large_val) = each %$hash;

    while (my ($key, $val) = each %$hash) {
        if ($val > $large_val) {
            $large_val = $val;
            $large_key = $key;
        }
    }
    return $large_key;
}

sub get_softclip_info{
    my $cigar = shift;
    my $seq   = shift;
    my $side  = shift;

    my $length     = 0;
    my $clipSeq    = '';

    if ( $side eq 'L' ){
	if ($cigar =~ m/^(\d+)S/){
	    $length = $1;
	}
	if($length > 0 ){
	    $clipSeq = substr($seq, 0, $length);
	}
    }
    elsif ( $side eq 'R'){
	if ($cigar =~ m/(\d+)S$/){
	    $length = $1;
	}
	if($length > 0 ){
	    $clipSeq = substr($seq, length($seq)-$length, $length);
	}
    }

    return($length, $clipSeq);

}


sub parse_strand_bitflag{
    my $bitflag = shift;
    my $strand  = '.';
    if($bitflag & 16){
        $strand  = '-';
    }
    else{
        $strand  = '+';
    }
    return($strand);
}

sub get_length_from_cigar {
  my $cigar_string = shift;
  my $cigar_length = 0;

  while($cigar_string=~m/(\d+)[MDX=N]/g){
    $cigar_length+=$1;
  }

  if($cigar_length == 0) {
    print STDERR "CIGAR string <'$cigar_string'> seems corrupt;\n";
  }
  return($cigar_length);
}

##############
## man page ##
##############

=pod

=head1 NAME

getMappingOverhang.pl [-l INT -c INT -cc INT] [--help|?] [--man] < FILE.sam

=head1 OPTIONS

=over 4

=item B<-l>

The minimal length of softclipped overhang to be considered

=item B<-c>

The minimal coverage at the overhang anchor site to be considered

=item B<-cc>

The minimal coverage within the overhang to be used for consensus sequence quality evaluation

=item B<--[help|?]>

Print the short version of the man-page.

=back

=head1 DESCRIPTION

Extracts soft clipped bases from an input sam file, constructs, characterizes and evaluates consensus sequence of the overhang.

Input sam file should single-end reads which are PCR artefact removed, unique, locally mapped. locally mapped, e.g., by STAR using

STAR --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --genomeDir ./STARIDX --readFilesIn reads.fa --alignEndsType Local --alignSoftClipAtReferenceEnds Yes

Output is a table, holding the following infos

=item B<chr> chromosome name

=item B<position>  chromosome position

=item B<strand>	strand [+-]

=item B<length>  length of the overhang

=item B<consensusSeq>  consensus sequence of the overhang

=item B<abs_conflictPosition>  number of positions in the overhang where at least one reads did not agree with the deduced consensus sequence

=item B<rel_conflictPosition>  ratio of positions in the overhang where at least one reads did not agree with the deduced consensus sequence

=item B<rel_BaseCountA>	ratio of bases in all overhang sequences which are A

=item B<rel_BaseCountConflicting>  ratio of bases in all overhang sequences which do not agree with the deduced consensus sequence (only positions with more than -cc INT reads are considered for Numerator and Denominator)

=head1 EXAMPLES

./getMappingOverhang.pl -l 4 -c 3 -cc 2 < FILE.sam > OUT.tsv

samtools view FILE.bam | ./getMappingOverhang.pl -l 4 -c 3 -cc 2 > OUT.tsv

=head1 AUTHOR

Fabian Amman, fabian@tbi.univie.ac.at

=cut
