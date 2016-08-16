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
                 'MAAFT' =>	1,
                 'RAXML' => 1,
                 'NOTUNG' => 1,
                 'NOISY' => 1,
                 'JAVA' => 1
               }

  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file,\%sections,scalar(@ARGV));
}




sub usage {
	return "USAGE: perl /path/to/run_genetrees.pl /path/to/config_file.ini";
}
