#!/usr/bin/perl -w

use strict;
use Cwd 'abs_path';
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
  'ENSEMBL' =>	{ 	'LOCAL' => 1
          },
  'DATABASE_CORE' =>	{ 	'HOST' => 1,
              'PORT' => 1,
              'RO_USER' => 1
            },
   'TAXA' =>	{},
   'SETUP' =>	{  'FASTA_DIR' => 1,
                 'REMOVE' =>	1,
                 'VIRTENV' => 1,
                 'ORTHOFINDER' => 1,
                 'BLASTPATH' => 1,
                 'MCL' => 1
               }

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


# read speciesID mapping from orthofinder directory if one exists
my %ids;
my $orthofinder_dir = $params->{'SETUP'}{'FASTA_DIR'}.'/Results/workingdirectory';
if (-d $orthofinder_dir){
  open IDS,"$orthofinder_dir/SpeciesIDs.txt";
  while (<IDS>){
    # add line to ids hash
  }
}

my %remove;
for my $taxon (@{$params->{'SETUP'}{'REMOVE'}}){
  # comment out matching proteomes and store ID
  system "perl -p -i -e 's/^[^#](.*$taxon)/#$1/'";
  $remove{$taxon} = $ids{$taxon};
}

for my $taxon (keys %{$params->{'TAXA'}}){
  # copy new proteomes to a new directory
  if (!$ids{$taxon} && !$remove{$taxon}){

  }
}

# prepare blastdbs with orthofinder

# run blast commands unless they refer to a removed proteome

# run orthofinder

# rename folder

sub usage {
	return "USAGE: perl /path/to/run_orthofinder.pl /path/to/config_file.ini";
}
