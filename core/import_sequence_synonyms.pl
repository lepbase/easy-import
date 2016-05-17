#!/usr/bin/perl -w

use strict;
use Cwd 'abs_path';
use File::Basename;

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
  'DATABASE_CORE' =>	{ 	'NAME' => 1,
              'HOST' => 1,
              'PORT' => 1,
              'RW_USER' => 1,
              'RW_PASS' => 1,
              'RO_USER' => 1
            },
  'META' =>	{	'ASSEMBLY.NAME' => 1},'
  FILES' => 	{	'SCAFFOLD_NAMES' => 1
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
foreach my $subsection (keys %{$params->{'FILES'}}){
	($infiles{$subsection}{'name'},$infiles{$subsection}{'type'}) = fetch_file($params->{'FILES'}{$subsection});
}

my $dbh = core_db_connect($params);


## load synonyms from file
## convert contig/scaffold names to standard format and load as synonyms
if ($params->{'FILES'}{'SCAFFOLD_NAMES'} || $params->{'SCAFFOLD_NAMES'}){
	add_seq_region_synonyms($params,\%infiles,$dbh);
}


sub usage {
	return "USAGE: perl /path/to/import_sequence_synonyms.pl /pat/to/config_file.ini";
}
