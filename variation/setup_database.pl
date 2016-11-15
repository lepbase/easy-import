#!/usr/bin/perl -w

use strict;
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
use EasyImport::Variation;

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

my $lib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl/modules';
my $variationlib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-variation/modules';
push @INC, $lib;
push @INC, $variationlib;

# create the variation database from a template

setup_variation_db($params);

sub usage {
	return "USAGE: perl setup_database.pl ini_file";
}
