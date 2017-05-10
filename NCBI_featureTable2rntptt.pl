#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;

my $featureTable_fh = "input_feature_table.txt.gz";
my $basename = '';

GetOptions(
  "ft=s"      => \$featureTable_fh,
  "man"       => sub{pod2usage(-verbose => 2)},
  "help|?"    => sub{pod2usage(-verbose => 1)}
    );

unless ( -e $featureTable_fh ){
  print STDERR "Info:\t$featureTable_fh was not found. Use parameter -ft to specify path to NCBI feature table\n";
  exit;
}

$basename = $featureTable_fh;
$basename =~s/_feature_table.txt.gz//;

my $rnt_fh          = "${basename}.rnt";
my $ptt_fh          = "${basename}.ptt";

open FT, "zcat $featureTable_fh | " or die "can t read $featureTable_fh via <zcat $featureTable_fh | .. >\n";
open RNT, "> $rnt_fh" or die "can t write to $rnt_fh\n";
open PTT, "> $ptt_fh" or die "can t write to $ptt_fh\n";

print STDERR "Processing $featureTable_fh\n";
print STDERR "Generating $rnt_fh and $ptt_fh\n";

print PTT "Species and strain info, complete genome. - 1..1234567\n";
print PTT "12345 proteins\n";
print PTT join("\t", qw/Location Strand Length PID Gene Synonym Code COG Product/)."\n";

print RNT "Species and strain info, complete genome. - 1..1234567\n";
print RNT "12345 RNAs\n";
print RNT join("\t", qw/Location Strand Length PID Gene Synonym Code COG Product/)."\n";


while(<FT>){
  chomp;
  next if (m/^#/);
  my @F = split"\t",$_;
  
  my $feature                 = ( defined($F[0]) && $F[0] )?($F[0]):("-");
  my $class                   = ( defined($F[1]) && $F[1] )?($F[1]):("-");
  my $assembly                = ( defined($F[2]) && $F[2] )?($F[2]):("-");
  my $assembly_unit           = ( defined($F[3]) && $F[3] )?($F[3]):("-");
  my $seq_type                = ( defined($F[4]) && $F[4] )?($F[4]):("-");
  my $chromosome              = ( defined($F[5]) && $F[5] )?($F[5]):("-");
  my $genomic_accession       = ( defined($F[6]) && $F[6] )?($F[6]):("-");
  my $start                   = ( defined($F[7]) && $F[7] )?($F[7]):("-");
  my $end                     = ( defined($F[8]) && $F[8] )?($F[8]):("-");
  my $strand                  = ( defined($F[9]) && $F[9] )?($F[9]):("-");
  my $product_accession       = ( defined($F[10]) && $F[10] )?($F[10]):("-");
  my $nonRedundant_refseq     = ( defined($F[11]) && $F[11] )?($F[11]):("-");
  my $related_accession       = ( defined($F[12]) && $F[12] )?($F[12]):("-");
  my $name                    = ( defined($F[13]) && $F[13] )?($F[13]):("-");
  my $symbol                  = ( defined($F[14]) && $F[14] )?($F[14]):("-");
  my $GeneID                  = ( defined($F[15]) && $F[15] )?($F[15]):("-");
  my $locus_tag               = ( defined($F[16]) && $F[16] )?($F[16]):("-");
  my $feature_interval_length = ( defined($F[17]) && $F[17] )?($F[17]):("-");
  my $product_length          = ( defined($F[18]) && $F[18] )?($F[18]):("-");
  my $attributes              = ( defined($F[19]) && $F[19] )?($F[19]):("-");


  if ($F[0] eq "CDS"){
    print PTT join("\t", "${start}..${end}", $strand, $product_length, $GeneID, $symbol, $locus_tag, "-", "-", $name)."\n";
  }
  elsif($F[0] eq "ncRNA" || $F[0] eq "rRNA" || $F[0] eq "tRNA" || $F[0] eq "tmRNA"){
    print RNT join("\t", "${start}..${end}", $strand, $feature_interval_length, $GeneID, $symbol, $locus_tag, "-", "-", $name)."\n";
  } 
}
close PTT;
close RNT;
close FT;

##############
## man page ##
##############

=pod

=head1 NAME

NCBI_featureTable2rntptt.pl -ft <FILE>

=head1 OPTIONS

=over 4

=item B<--ft> <FILE>

NCBI feature table gzipped.

=item B<--[help|?]>

Print the short version of the man-page.

=item B<--[man]>

Print man-page.

=back

=head1 DESCRIPTION

Usefull to produce old style *.rnt *.ptt files from newer file format feature table as currently provided by NCBI.

=head1 EXAMPLES

perl NCBI_featureTable2rntptt.pl -ft GCF_000005845.2_ASM584v2_feature_table.txt.gz

=head1 AUTHOR

Fabian Amman, fabian@tbi.univie.ac.at

=cut
