#!/usr/bin/perl -w

use strict;
use strict;
use Cwd 'abs_path';
use File::Basename;
use File::Find;
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
  'ENSEMBL' =>	{ 	'LOCAL' => 1
          },
  'DATABASE_COMPARA' =>	{ 	'NAME' => 1,
              'HOST' => 1,
              'PORT' => 1,
              'RW_USER' => 1,
              'RW_PASS' => 1,
              'RO_USER' => 1
            },
    'TAXA' => {},
    'ORTHOGROUP' => {  'PREFIX' => 1,
  	           'PROTEIN' => 1,
  	           'PROTEIN_ALIGN' => 1,
  	           'PROTEIN_TRIMMED' => 1,
  	           'FNAFILE' => 1,
  	           'BOUNDEDFILE' => 1,
  	           'TREE' => 1
            }
  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file,\%sections);
}

my $lib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl/modules';
my $comparalib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-compara/modules';
push @INC, $lib;
push @INC, $comparalib;
load Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
load Bio::EnsEMBL::Compara::Graph::NewickParser;
load Bio::EnsEMBL::Compara::GeneTreeNode;
load Bio::EnsEMBL::Compara::GeneTreeMember;
load Bio::EnsEMBL::Compara::SpeciesTreeNode;
load Bio::EnsEMBL::Compara::SpeciesTree;
load Bio::EnsEMBL::Compara::DBSQL::NCBITaxonAdaptor;

my $dbh = compara_db_connect($params);

our $prefix = $params->{'ORTHOGROUP'}{'PREFIX'};

my $taxlist;
foreach my $key (keys %{$params->{'TAXA'}}){
  $taxlist .= $key.'|';
}
chop $taxlist;
$params->{'ORTHOGROUP'}{'TAXA'} = $taxlist;

my ($st_nodes,$core_dbs) = fetch_species_tree_nodes($params,$dbh);

find({wanted => sub {
  my $file = $File::Find::name;
  if (!-d $file){
    if ($file =~ m/$prefix\w+$/){
      warn "importing $file\n";
      load_sequences($dbh,$params,$st_nodes,$core_dbs,$file);
    }
  }
},
no_chdir => 1},'.');

sub usage {
        return "USAGE: perl -I /path/to/dir/containing/Ensembl_Import.pm /path/to/import_homologues.pl ini_file";
}
