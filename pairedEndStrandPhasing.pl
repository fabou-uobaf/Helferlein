#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;

## Parameter
my $bams_fh = '';	# holds ',' seperated list of inout bam files
my @bams_fh = ();	# array holding input bam files
my @out_fh  = ();	# generated output file names 
my $samtools = "";	# path to samtools
my $switchfirst = 0;    # flag if strand if first in pair should be switched
my $switchsecond = 0;	# flag if strand if second in pair should be switched
my $switchSE = 0;       # flag if strand if single end aln should be switched
my $ppFilter = 0;       # flag if only proper paired reads should be considered
my %actionCounter = (); # log number of alignments which are changed/not changed

GetOptions(
  "bams=s"        => \$bams_fh,
  "samtools=s"    => \$samtools,
  "switchfirst"   => \$switchfirst,
  "switchsecond"  => \$switchsecond,
  "switchSE"      => \$switchSE,
  "man"           => sub{pod2usage(-verbose => 2)},
  "help|?"        => sub{pod2usage(-verbose => 1)}
    );

# check input variables
@bams_fh = split",", $bams_fh;
if(@bams_fh){
  foreach my $fh (@bams_fh){
    if (-e $fh){
      my $out_fh = $fh;
      $out_fh=~s/\.bam/.strandPhased.bam/;
      push @out_fh, $out_fh;
      print STDERR "Info:\t$fh was found. Strand phased output will be written to $out_fh\n";
    }
    else{
      print STDERR "Error:\tfile $fh does not exist\t=> aborted!\n";
      exit;
    }
  }
}
else{
  print STDERR "Error:\tno input files specified. Please use option --bams <FILE1,FILE2,..>\n";
}

if ($samtools){
  if(-e $samtools){
    print STDERR "Info:\tsamtools at $samtools will be used\n";
  }
  elsif (my $bin = `which samtools`){
    chomp($bin);
    $samtools = $bin;
    print STDERR "Info:\tsamtools at $samtools will be used\n";
  }
  else{
    print STDERR "Error:\tno samtools found. Please use option --samtools <PATH>\n";
  }
}
else{
  if (my $bin = `which samtools`){
    chomp($bin);
    $samtools = $bin;
    print STDERR "Info:\tsamtools at $samtools will be used\n"; 
  }
  else{
    print STDERR "Error:\tno samtools found. Please use option --samtools <PATH>\n";
  }
}

if ($switchfirst || $switchsecond || $switchSE){
  print STDERR "Info:\tStand of first in read alignments will be switched.\n" if ($switchfirst);
  print STDERR "Info:\tStand of second in read alignments will be switched.\n" if ($switchsecond);
  print STDERR "Info:\tStand of single-end read alignments will be switched.\n" if ($switchSE);
}
else{
  print STDERR "Info:\tStand of no alignments will be switched. Use parameters --switchfirst, --switchsecond or --switchSE to specify reads to be strand switched\n";
}

# process input bam files
foreach my $idx (0..$#bams_fh){
  my $fh = $bams_fh[$idx];
  my $out = $out_fh[$idx];

  my $cmdO = "| $samtools view -bS - | $samtools sort - -T foo > $out";
  open O, $cmdO or die "Error:\tcan't write bam file via <$cmdO>\n";
  my $cmdI = "$samtools view -h $fh | ";
  open I, $cmdI or die "Error:\tcan't read bam file via <$cmdI>\n";
  while(<I>){
    chomp;
    if(m/^@/){
      print O "$_\n";
      next;
    }
    else{
      my @F=split"\t", $_;
      my $bitFlag = $F[1];
      my $FinR = ($bitFlag & 64)?(1):(0);
      my $SinR = ($bitFlag & 128)?(1):(0);
      my $SE   = ($bitFlag & 1)?(0):(1);

      if ($switchfirst && $FinR){
        if ($bitFlag & 16){
          $bitFlag -= 16;
          $actionCounter{$out}->{firstinread}->{rev2fwd}++ if ($FinR);
        }
        else{
          $bitFlag += 16;
          $actionCounter{$out}->{firstinread}->{fwd2rev}++ if ($FinR);
        }
      }
      else{
        $actionCounter{$out}->{firstinread}->{unchanged}++ if ($FinR);
      }

      if ($switchsecond && $SinR){
        if ($bitFlag & 16){
          $bitFlag -= 16;
          $actionCounter{$out}->{secondinread}->{rev2fwd}++ if ($SinR); 
        }   
        else{
          $bitFlag += 16;
          $actionCounter{$out}->{secondinread}->{fwd2rev}++ if ($SinR);
        }   
      }
      else{
        $actionCounter{$out}->{secondinread}->{unchanged}++ if ($SinR);
      }

      if ($switchSE && $SE){
        if ($bitFlag & 16){
          $bitFlag -= 16;
          $actionCounter{$out}->{singleendread}->{rev2fwd}++ if ($SE);
        }
        else{
          $bitFlag += 16;
          $actionCounter{$out}->{singleendread}->{fwd2rev}++ if ($SE);
        }
      }
      else{
        $actionCounter{$out}->{singleendread}->{unchanged}++ if ($SE);
      }


      $F[1] = $bitFlag;
      print O join("\t", @F)."\n";
    }
  }
  close I;
  close O;
}

foreach my $sample (keys %actionCounter){
  print STDERR "Stats:\t$sample\n";
  foreach my $type (qw/firstinread secondinread singleendread/){
    print STDERR "\t\t\t\t$type\n";
    foreach my $action (qw/unchanged fwd2rev rev2fwd/){
      my $count = 0;
      if ( defined($actionCounter{$sample}->{$type}->{$action}) ){
        $count = $actionCounter{$sample}->{$type}->{$action};
      }
      my $aaction = $action;
      $aaction .="\t" unless ( $aaction eq "unchanged");
      print STDERR "\t\t\t\t\t$aaction\t$count\n";
    }
  }
}

##############
## man page ##
##############

=pod

=head1 NAME

pairedEndStrandPhasing.pl --bams <BAM1>,<BAM2>,...,<BAMn> [--switchfirst] [--switchsecond] [--switchSE] [--samtools <PATH>] [--man] [--help]

=head1 OPTIONS

=over 4

=item B<--bams> <FILE>

Comma-separated list of bam files to be processed.

=item B<--switchfirst> 

Flag if alignments classified as first in pair should be strand switched (Default: no action)

=item B<--switchsecond>

Flag if alignments classified as second in pair should be strand switched (Default: no action)

=item B<--switchSE>

Flag if alignments classified as single-end reads (no read pair) should be strand switched (Default: no action)

=item B<--samtools> <STRING>

Path to samtools which are used to read and write bam files. Tested with samtools version 1.3 and 1.4. (Default: samtools specified in the \$PATH environment.

=item B<--[help|?]>

Print the short version of the man-page.

=item B<--[man]>

Print man-page.


=back

=head1 DESCRIPTION

Switches the strand of alignments in a bam file. Can treat pair/mates and single-end reads differently. 

Useful if sequencing libraries are of type "1+-,1-+,2++,2--" (Lexogene library) or "1++,1--,2+-,2-+" (default Illumina library) and follow up tools consider pair/mates as single reads with their respective strand information.

=head1 EXAMPLES

perl ./pairedEndStrandPhasing.pl --bams Ecoli_rep1.bam,Ecoli_rep2.bam --switchfirst --switchsecond --switchSE --samtools ${HOME}/.local/bin/samtools

=head1 AUTHOR

Fabian Amman, fabian@tbi.univie.ac.at

=cut
