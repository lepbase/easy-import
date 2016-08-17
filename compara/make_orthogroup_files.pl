#!/usr/bin/env perl

use strict;
use warnings;

my $orthologousgroupstxt_filename = shift @ARGV;

my @species_list = @ARGV;
my $species_regexp = join("|", @species_list);
# or get this from params ini file like import_orthogroups.pl

my %sequences;

# load all sequences in memory
for my $species (@species_list) {
  $sequences{$species}{faa} = fastafile2hash("$species\_-_canonical_proteins.fa");
  $sequences{$species}{fna} = fastafile2hash("$species\_-_cds_translationid.fa");
  $sequences{$species}{fba} = fastafile2hash("$species\_-_protein_bounded_exon.fa");
}

# create folder with sequences for each orthogroup
open OG, "<$orthologousgroupstxt_filename" or die $!;
while (<OG>) {
  my @tokens = split /\s+/;
  my $orthogroup_id = shift @tokens;
  $orthogroup_id =~ s/:$//;
  mkdir "orthogroups/$orthogroup_id";
  open IDS, ">orthogroups/$orthogroup_id/$orthogroup_id" or die $!;
  open FAA, ">orthogroups/$orthogroup_id/$orthogroup_id.faa" or die $!;
  open FNA, ">orthogroups/$orthogroup_id/$orthogroup_id.fna" or die $!;
  open FBA, ">orthogroups/$orthogroup_id/$orthogroup_id.fba" or die $!;
  for my $sequence_id (@tokens) {
    if ($sequence_id =~ /^(${species_regexp})_\S+$/) {
      my $species = $1;
      if (exists $sequences{$species}{faa}{$sequence_id} and 
          exists $sequences{$species}{fna}{$sequence_id} and
          exists $sequences{$species}{fba}{$sequence_id}) {
        print IDS "$sequence_id\n";
        print FAA ">$sequence_id\n" . $sequences{$species}{faa}{$sequence_id}->{seq} . "\n";
        print FNA ">$sequence_id\n" . $sequences{$species}{fna}{$sequence_id}->{seq} . "\n";
        print FBA ">$sequence_id\n" . $sequences{$species}{fba}{$sequence_id}->{seq} . "\n";
      }
      else {
        warn "$sequence_id does not exist in one of the input fasta files";
      }
    }
    else {
      warn "$sequence_id does not match species list";
    }
  }
  close IDS;
  close FAA;
  close FNA;
  close FBA;
}

#############################################################################
sub fastafile2hash
{
  my $fastafile  = shift @_;
  my $changecase = "N";
  my $order      = "S"; # S = same as input, or R = random
  $changecase    = substr(uc(shift @_),0,1) if @_;
  $order         = substr(uc(shift @_),0,1) if @_;
  my %sequences;
  my $fh = &read_fh($fastafile);
  my $seqid;
  my $seq_counter;
  while (<$fh>)
  {
    if (/^>(\S+)(.*)/) {
      $seqid = $1;
      $sequences{$seqid}{desc} = $2;
      $sequences{$seqid}{order} = $order eq "S" ? $seq_counter++ : rand;
    }
    else {
      if (/\d/) {
        chomp($sequences{$seqid}{seq} .= " $_"); # add space to sep qual values
        $sequences{$seqid}{seq} =~ s/^\s+//;
        $sequences{$seqid}{seq} =~ s/\s+$//;
        next;
      }
      chomp($sequences{$seqid}{seq} .= lc($_)) if $changecase eq "L";
      chomp($sequences{$seqid}{seq} .= uc($_)) if $changecase eq "U";
      chomp($sequences{$seqid}{seq} .= $_    ) if $changecase eq "N";
    }
  }
  return \%sequences;
}

#############################################################################

sub read_fh {
  my $filename = shift @_;
  my $filehandle;
  if ($filename =~ /gz$/) {
    open $filehandle, "gunzip -dc $filename |" or die $!;
  }
  else {
    open $filehandle, "<$filename" or die $!;
  }
  return $filehandle;
}

