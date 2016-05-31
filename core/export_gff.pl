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
  'DATABASE_CORE' =>	{ 	'NAME' => 1,
              'HOST' => 1,
              'PORT' => 1,
              'RW_USER' => 1,
              'RW_PASS' => 1,
              'RO_USER' => 1
            }
  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file,\%sections);
}

## download/obtain files using methods suggested by file paths and extensions
my %infiles;
foreach my $subsection (sort keys %{$params->{'FILES'}}){
	($infiles{$subsection}{'name'},$infiles{$subsection}{'type'}) = fetch_file($params->{'FILES'}{$subsection});
}

my $lib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl/modules';
my $iolib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-io/modules';
my $comparalib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-compara/modules';
push @INC, $lib;
push @INC, $iolib;
push @INC, $comparalib;
load Bio::EnsEMBL::Registry;
load Bio::EnsEMBL::Utils::IO::GFFSerializer;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => $params->{'DATABASE_CORE'}{'HOST'},
    -user => $params->{'DATABASE_CORE'}{'RO_USER'},
    -port   => $params->{'DATABASE_CORE'}{'PORT'},
    -driver => 'mysql',
);

my $dbname = $params->{'DATABASE_CORE'}{'NAME'};

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => $params->{'DATABASE_CORE'}{'RO_USER'},
    -dbname => $dbname,
    -host   => $params->{'DATABASE_CORE'}{'HOST'},
    -port   => $params->{'DATABASE_CORE'}{'PORT'},
    -driver => 'mysql'
);


my $meta_container = $dba->get_adaptor("MetaContainer");

my $display_name    = $meta_container->get_display_name();

my   $output_fh;
open $output_fh, ">$display_name.gff";

my   $ontology_adaptor = $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
my   $serializer       = Bio::EnsEMBL::Utils::IO::GFFSerializer->new($ontology_adaptor,$output_fh);

my   $slice_adaptor    = $dba->get_adaptor("Slice");
my   $supercontigs     = $slice_adaptor->fetch_all('toplevel');

$serializer->print_main_header($supercontigs, $slice_adaptor, $output_fh);

foreach my $slice (@{$supercontigs}) {
    $serializer->print_feature_list($slice->get_all_Genes);
}


sub usage {
	return "USAGE: perl /path/to/export_gff.pl /path/to/config_file.ini";
}
