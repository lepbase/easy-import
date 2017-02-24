#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';
use File::Path qw(make_path);
use File::Basename;
use Module::Load;
use File::Find;

## find the full path to the directory that this script is executing in
our $dirname;
BEGIN {
  $dirname  = dirname(abs_path($0));
}
use lib "$dirname/../modules";
use lib "$dirname/../gff-parser";
use EasyImport::Core;
use EasyImport::Compara;

## load parameters from an INI-style config file
my %sections = (
  'ENSEMBL' =>  {
    'LOCAL' => 1,
  },
  'TAXA' => {},
  'SETUP' => {
    'FASTA_DIR' => 1,
    'MAFFT' => 1,
    'NOISY' => 1,
    'RAXML' => 1,
    'NOTUNG' => 1,
    'NOTUNG_SPECIESTREE' => 1,
    'JAVA' => 1,
  },
  'ORTHOGROUP' => {
    'PREFIX' => 1,
    'SUFFIXLENGTH' => 1,
    'PROTEIN' => 1,
    'PROTEIN_ALIGN' => 1,
    'PROTEIN_TRIMMED' => 1,
    'FNAFILE' => 1,
    'BOUNDEDFILE' => 1,
    'TREE' => 1,
    'HOMOLOG' => 1,
  },
  'SPECIES_SET' => {
    'TREE_FILE' => 1,
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
my $comparalib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-compara/modules';
push @INC, $lib;
push @INC, $comparalib;

#===============================================================================

my $orthogroup_prefix       = $params->{'ORTHOGROUP'}{'PREFIX'};

open QSUB, ">run_genetrees_qsub.txt" or die $!;

find({wanted => sub {
  my $file = $File::Find::name;
  if (!-d $file){
    if ($file =~ m/$orthogroup_prefix\w+$/){
      warn "creating genetree qsub script for $file\n";
      create_qsub_script($params, $file);
      print QSUB "$file.bash\n";
    }
  }
},
no_chdir => 1},'.');

close QSUB;
open  QSUB, ">qsub.bash" or die $!;
print QSUB  '
#$ -j y
#$ -pe smp 2
export MAFFT='. $params->{'SETUP'}{'MAFFT'}.'
export RAXML='. $params->{'SETUP'}{'RAXML'}.'
export NOISY='. $params->{'SETUP'}{'NOISY'}.'
export NOTUNG='. $params->{'SETUP'}{'NOTUNG'}.'
export NOTUNG_SPECIESTREE='. $params->{'SETUP'}{'NOTUNG_SPECIESTREE'}.'
export COMPARATEMP='. $params->{'SETUP'}{'COMPARATEMP'}.'
export JAVA='. $params->{'SETUP'}{'JAVA'}.'
awk "NR==$SGE_TASK_ID" run_genetrees_qsub.txt | bash
';

#system "qsub -l h=bigshot|bigwig|bigbird|bigbang -cwd -V -t 1-`cat run_genetrees_qsub.txt | wc -l` qsub.bash";

sub usage {
  return "USAGE: perl /path/to/run_genetrees.pl /path/to/config_file.ini";
}
