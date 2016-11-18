#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';
use File::Path qw(make_path);
use File::Basename;
use Module::Load;

## find the full path to the directory that this script is executing in
our $dirname;
BEGIN {
  $dirname  = dirname(abs_path($0));
}
use lib "$dirname/../modules";
use lib "$dirname/../gff-parser";
use EasyImport::Core;

## load parameters from an INI-style config file
my %sections = (
  'ENSEMBL' =>  {
    'LOCAL' => 1,
  },
  'TAXA' =>  {},
  'SETUP' => {
    'FASTA_DIR' => 1,
    'REMOVE'    =>  1,
  },
  'ORTHOGROUP'  => {
    'PREFIX'       => 1,
    'SUFFIXLENGTH' => 1,
  },
);

## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
  load_ini($params,$ini_file,\%sections,scalar(@ARGV));
}

my $lib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl/modules';
push @INC, $lib;
load Bio::EnsEMBL::DBSQL::DBAdaptor;

my $fastadir = $params->{'SETUP'}{'FASTA_DIR'};
die "Core sequences fasta dir $fastadir does not exist\n" unless -d $fastadir;

my $orthologousgroupstxt_filename = $params->{'ORTHOGROUP'}{'ORTHOGROUPS_FILE'};
die "Orthogroups file $orthologousgroupstxt_filename does not exist\n" unless -s $orthologousgroupstxt_filename;

my $orthogroup_prefix       = $params->{'ORTHOGROUP'}{'PREFIX'};
my $orthogroup_suffixlength = $params->{'ORTHOGROUP'}{'SUFFIXLENGTH'};

my @species_list = keys %{$params->{'TAXA'}};
my $species_regexp = join("|", @species_list);

my %sequences;

# load all sequences in memory
for my $species (@species_list) {
  $sequences{$species}{faa} = fastafile2hash("$fastadir/canonical_proteins/$species\_-_canonical_proteins.fa");
  $sequences{$species}{fna} = fastafile2hash("$fastadir/canonical_cds_translationid/$species\_-_canonical_cds_translationid.fa");
  $sequences{$species}{fba} = fastafile2hash("$fastadir/canonical_protein_bounded_exon/$species\_-_canonical_protein_bounded_exon.fa");
}

# create folder with sequences for each orthogroup
open OG, "<$orthologousgroupstxt_filename" or die $!;
while (<OG>) {
  my @tokens = split /\s+/;
  my $orthogroup_id = shift @tokens;
  next if scalar @tokens < 4; # no need to make trees if there are less than 4 sequences
  $orthogroup_id =~ s/.*?(.{$orthogroup_suffixlength}):$/$1/;
  $orthogroup_id = $orthogroup_prefix . $orthogroup_id;
  # In the future, add checks to ensure that an orthogroup ID stays consistent across releases
  make_path "orthogroups/$orthogroup_id";
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

