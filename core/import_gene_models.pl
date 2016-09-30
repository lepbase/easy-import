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

## load parameters from an INI-style config file for maximal ensemblism
my %sections = (
  'DATABASE_CORE' =>	{ 	'NAME' => 1,
              'HOST' => 1,
              'PORT' => 1,
              'RW_USER' => 1,
              'RW_PASS' => 1,
              'RO_USER' => 1
            },
  'FILES' => 	{	'GFF' => 1
        }
  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file,\%sections,scalar(@ARGV));
}

## connect to core database
my $dbh = core_db_connect($params);

## download/obtain files using methods suggested by file paths and extensions
my %infiles;
foreach my $subsection (sort keys %{$params->{'FILES'}}){
	($infiles{$subsection}{'name'},$infiles{$subsection}{'type'}) = fetch_file($params->{'FILES'}{$subsection});
}

## truncate database tables if option is specified
if ($params->{'MODIFY'}{'TRUNCATE_GENE_TABLES'}){
	truncate_gene_tables($dbh);
}

## load gff into ensembl
my $suffix = '.gff';
$suffix = '.sorted'.$suffix if $params->{'GFF'}{'SORT'};
my $filename = $infiles{'GFF'}{'name'}.$suffix;
if (!-e $filename){
  $filename = $infiles{'GFF'}{'name'};
  warn "WARNING: unable to locate prepared GFF, attempting to use unmodified file\n"
}

gff_to_ensembl($filename,$dbh,$params);

count_rows($dbh,qw( gene transcript translation exon));


sub usage {
	return "USAGE: perl /path/to/import_gene_models.pl /pat/to/config_file.ini";
}
