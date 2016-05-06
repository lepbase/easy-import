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
  'ENSEMBL' =>	{ 	'LOCAL' => 1
          },
  'DATABASE_CORE' =>	{ 	'NAME' => 1,
              'HOST' => 1,
              'PORT' => 1,
              'RW_USER' => 1,
              'RW_PASS' => 1,
              'RO_USER' => 1
            },
  'DATABASE_TAXONOMY' => {	'NAME' => 1,
                'HOST' => 1,
                'PORT' => 1,
                'RO_USER' => 1
              },
  'DATABASE_TEMPLATE' => 	{	'NAME' => 1,
                'HOST' => 1,
                'PORT' => 1,
                'RO_USER' => 1
              },
  'META' =>	{	'SPECIES.PRODUCTION_NAME' => 1,
          'SPECIES.SCIENTIFIC_NAME' => 1,
          'SPECIES.TAXONOMY_ID' => 1,
          'ASSEMBLY.NAME' => 1,
          'GENEBUILD.METHOD' => 1,
          'PROVIDER.NAME' => 1,
          'PROVIDER.URL' => 1
        },
  'FILES' => 	{	'SCAFFOLD' => [ 'CONTIG' ],
          'CONTIG' => [ 'SCAFFOLD']
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

if (table_exists( $dbh, "meta") && !$params->{'MODIFY'}{'OVERWRITE_DB'}){
	if ($params->{'MODIFY'}{'TRUNCATE_SEQUENCE_TABLES'}){
		truncate_seq_tables($dbh);
	}
	else {
		## print a warning message saying that running this script may not have been a good idea
		warn "WARNING: Database exists and [MODIFY] TRUNCATE_TABLES was not set, unexpected things may happen.\n";
		warn "WARNING: If this was not what you intended you have 5 seconds to hit CTRL-C...\n";
		sleep(5);
	}
}
else {
	## create database, load schema and populate meta tables etc.
	setup_core_db($dbh,$params);
}

## load contigs/scaffolds and agp
## TODO: handle levels other than scaffold and contig
load_sequences($params,\%infiles,$dbh);


## load synonyms from file
## convert contig/scaffold names to standard format and load as synonyms
if ($params->{'FILES'}{'SCAFFOLD_NAMES'} || $params->{'SCAFFOLD_NAMES'}){
	add_seq_region_synonyms($params,\%infiles,$dbh);
}

## count the number of sequences that have been loaded into the db
## TODO: check that this is the expected number of sequences
count_sequences($params,\%infiles,$dbh);


sub usage {
	return "USAGE: perl /path/to/import_sequences.pl /pat/to/config_file.ini";
}
