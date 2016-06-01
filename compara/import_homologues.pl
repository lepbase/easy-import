#!/usr/bin/perl -w

use strict;
use strict;
use Cwd 'abs_path';
use File::Basename;
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
  'ENSEMBL' =>	{ 	'LOCAL' => 1
          },
  'DATABASE_COMPARA' =>	{ 	'NAME' => 1,
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

my $lib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl/modules';
my $comparalib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-compara/modules';
push @INC, $lib;
push @INC, $comparalib;

my $dbh = compara_db_connect($params);

opendir DIR, $params->{'ORTHOGROUP'}{'PATH'} or die "Cannot open directory: $!";
my $prefix = $params->{'ORTHOGROUP'}{'PREFIX'};
my @files = grep { /^$prefix\d+$/ } readdir DIR;
closedir DIR;

find(&wanted,$params->{'ORTHOGRUP'}{'PATH'})

sub wanted {
  my $file = shift;
  my $prefix = $params->{'ORTHOGRUP'}{'PREFIX'};
  if ($file =~ m/^$prefix\w+$/){
    load_sequences($dbh,$params,$file);
  }
}


sub usage {
	return "USAGE: perl -I /path/to/dir/containing/Ensembl_Import.pm /path/to/import_homologues.pl ini_file";
}
