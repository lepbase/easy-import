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
use Ensembl_Import;

## load parameters from an INI-style config file for maximal ensemblism
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file);
}

## connect to core database
my $dbh = core_db_connect($params);

## download/obtain files using methods suggested by file paths and extensions
my %infiles;
foreach my $subsection (keys %{$params->{'FILES'}}){
	($infiles{$subsection}{'name'},$infiles{$subsection}{'type'}) = fetch_file($params->{'FILES'}{$subsection});
}

## truncate database tables if option is specified
if ($params->{'MODIFY'}{'TRUNCATE_GENE_TABLES'}){
	truncate_gene_tables($dbh);
}

## load gff into ensembl
my $suffix = '.gff';
$suffix = '.sorted'.$suffix if $params->{'GFF'}{'SORT'};

gff_to_ensembl($infiles{'GFF'}{'name'}.$suffix,$dbh,$params);

count_rows($dbh,qw( gene transcript translation exon));


sub usage {
	return "USAGE: perl -I /path/to/dir/containing/Ensembl_Import.pm /path/to/gene_model_import.pl ini_file";
}
