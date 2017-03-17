#!/usr/bin/perl

use strict;
use warnings;


my ($removed, $kept, $total) = (0,0,0);
my ($prevChr, $prevPos, $prevStrand, $prevSeq, $prevCIGAR) = (-1,-1,-1,-1,-1);
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
  if ($prevChr eq $currChr && $prevPos eq $currPos && $prevStrand eq $currStrand && $prevSeq eq $currSeq && $prevCIGAR eq $currCIGAR){
    $removed++;
  }
  else{
    print join("\t", @F)."\n";
    $kept++;
    ($prevChr, $prevPos, $prevStrand, $prevSeq, $prevCIGAR) = ($currChr, $currPos, $currStrand, $currSeq, $currCIGAR);
  }
}

if ($total){
  printf STDERR "Total reads:\t%d\n", $total;
  printf STDERR "Reads kept:\t%d\t(%.3f)\n", $kept, $kept/$total;
  printf STDERR "Reads removed:\t%d\t(%.3f)\n", $removed, $removed/$total;
}
else{
  print STDERR "No reads in the input file\n";
}
