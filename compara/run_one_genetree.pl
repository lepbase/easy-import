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
use EasyImport::Compara;

## load parameters from an INI-style config file
my %sections = (
  'ENSEMBL' =>	{  'LOCAL' => 1
                },
  'DATABASE_CORE' =>	{  'HOST' => 1,
                         'PORT' => 1,
                         'RO_USER' => 1
                      },
   'TAXA' =>	{},
   'SETUP' =>	{  'FASTA_DIR' => 1,
                 'MAFFT' =>	1,
                 'NOISY' => 1,
                 'RAXML' => 1,
                 'NOTUNG' => 1,
                 'JAVA' => 1
               }

  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one orthogroup eg OG0012345\n",usage(),"\n" unless $ARGV[0];
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[1];
my $orthogroup_id = shift @ARGV;
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
  load_ini($params,$ini_file,\%sections,scalar(@ARGV));
}
my $orthogroup_prefix;
die "ini files seem incorrect\n" unless $orthogroup_prefix = $params->{'ORTHOGROUP'}{'PREFIX'};
die "$orthogroup_id does not match pattern " . $params->{'ORTHOGROUP'}{'PREFIX'} . '\w+' unless $orthogroup_id =~ /$orthogroup_prefix\w+/;

my $lib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl/modules';
my $comparalib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-compara/modules';
push @INC, $lib;
push @INC, $comparalib;

print "cd orthogroups/$orthogroup_id;";
if (exists $params->{'SETUP'}{'MAFFT'}) {
  print $params->{'SETUP'}{'MAFFT'} . " --treeout --auto --reorder $orthogroup_id.faa > $orthogroup_id.faa.mafft;" ; 
}
if (exists $params->{'SETUP'}{'MAFFT'} and
    exists $params->{'SETUP'}{'NOISY'}) {
  print $params->{'SETUP'}{'NOISY'} . "  --seqtype P $orthogroup_id.faa.mafft;" ; 
}
if (exists $params->{'SETUP'}{'MAFFT'} and
    exists $params->{'SETUP'}{'NOISY'} and
    exists $params->{'SETUP'}{'RAXML'}) {
  print $params->{'SETUP'}{'RAXML'} . " -f a -x 12345 -# 3 -T 1 -p 12345 -m PROTGAMMAAUTO -s $orthogroup_id.faa_out.fas -n $orthogroup_id;";
  print "mv RAxML_info.$orthogroup_id $orthogroup_id.RAxML_info;";
  print "mv RAxML_bipartitionsBranchLabels.$orthogroup_id $orthogroup_id.RAxML_bipartitionsBranchLabels;";
  print "mv RAxML_bipartitions.$orthogroup_id $orthogroup_id.RAxML_bipartitions;";
  print "mv RAxML_bestTree.$orthogroup_id $orthogroup_id.RAxML_bestTree;";
  print "mv RAxML_bootstrap.$orthogroup_id $orthogroup_id.RAxML_bootstrap;";
}
if (exists $params->{'SETUP'}{'MAFFT'} and
    exists $params->{'SETUP'}{'NOISY'} and
    exists $params->{'SETUP'}{'RAXML'} and
    exists $params->{'SETUP'}{'NOTUNG'} and
    exists $params->{'SETUP'}{'NOTUNG_SPECIESTREE'}) {
  print "java -jar " . $params->{'SETUP'}{'NOTUNG'} . 
    " --treeoutput nhx --root -s " . 
    $params->{'SETUP'}{'NOTUNG_SPECIESTREE'} . 
    " -g $orthogroup_id.RAxML_bipartitionsBranchLabels;";
  print "java -jar " . $params->{'SETUP'}{'NOTUNG'} . 
    " --nolosses --treeoutput nhx --homologtabletabs --reconcile --stpruned -s " . 
    $params->{'SETUP'}{'NOTUNG_SPECIESTREE'} . 
    " -g $orthogroup_id.RAxML_bipartitionsBranchLabels.rooting.0;";
}

sub usage {
  return "USAGE: perl /path/to/run_one_genetree.pl OG0012345 /path/to/config_file.ini";
}
