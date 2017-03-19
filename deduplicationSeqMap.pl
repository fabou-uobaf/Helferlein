#!/usr/bin/perl

use strict;
use warnings;


my ($removed, $kept, $samePos, $total) = (0,0,0,0);
my ($prevChr, $prevPos, $prevStrand, $prevSeq, $prevCIGAR) = (-1,-1,-1,-1,-1);
my %entry = ();

while(<>){
  chomp;
 
  # print header
  if(m/^@/){
    print "$_\n";
    next;
  }
  
  my @F=split"\t", $_;
  $total++;

  my ($currChr, $currPos, $currStrand, $currSeq, $currCIGAR) = ($F[2],$F[3],$F[1],$F[9],$F[5]);
  
  # check if input is sorted
  if ($currStrand eq $prevStrand && $currPos < $prevPos){
    print STDERR "Aborted:\tThe input sam files seems not be positionally sorted\n";
    exit;
  }
  
  # check if same mapping position of current and previous read
  if ($prevChr eq $currChr && $prevPos == $currPos) {
    $entry{$currStrand}->{$currSeq}->{$currCIGAR} = $_;
  }
  else{
    $samePos++;
    $kept += &report(\%entry);
    %entry=();
    $entry{$currStrand}->{$currSeq}->{$currCIGAR} = $_;
    ($prevChr, $prevPos, $prevStrand, $prevSeq, $prevCIGAR) = ($currChr, $currPos, $currStrand, $currSeq, $currCIGAR);
  }
}
$samePos++;
$kept += &report(\%entry);

if ($total){
  $removed = $total - $kept;
  printf STDERR "Total alignments:\t%d\n", $total;
  printf STDERR "Alignments kept:\t%d\t(%.3f)\n", $kept, $kept/$total;
  printf STDERR "Alignments removed:\t%d\t(%.3f)\n", $removed, $removed/$total;
  printf STDERR "Uniq alignment starts:\t%d\t(%.3f)\n", $samePos, $samePos/$total;
}
else{
  print STDERR "No reads in the input file\n";
}

## subroutines
sub report {
  my $h = shift;
  my $n = 0;
  
  foreach my $k1 (sort keys %$h){
    foreach my $k2 (sort keys %{$h->{$k1}}){
      foreach my $k3 (sort keys %{$h->{$k1}->{$k2}}){
        print "$h->{$k1}->{$k2}->{$k3}\n";
        $n++;
      } 
    }
  }
  return($n);
}
