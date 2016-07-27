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
use EasyImport::Search;

# create indexes for autocomplete and for search, both on descriptions and on single word terms
# for autocomplete, order will be
#		single word LIKE with wildcard at end
#		description MATCH with full text index
#		description LIKE with wildcard at start and end
# results will be returned as soon as there are at least 10 matches
# strings with and without spaces will be duplicated in a separate tables to autocomplete faster
# table structure will be
#		multi - all search terms as varchar(255)
#		multi_32 - single word terms as varchar(32)
#		multi_255 - multi word terms as varchar(255)
# each table will be duplicated for single production name searches to make single species search faster
#		production_name - all search terms as varchar(255)
#		production_name_32 - single word terms as varchar(32)
#		production_name_255 - multi word terms as varchar(255)


## load parameters from an INI-style config file
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file);
}

## connect to the core_db
my $core_dbh = core_db_connect($params);

# connect to/create search database
my $search_dbh = search_db_connect($params,$dirname);

my $production_name = lc $params->{'META'}{'SPECIES.PRODUCTION_NAME'};

# remove any existing data
$search_dbh->do("DROP TABLE IF EXISTS $production_name");
$search_dbh->do("DROP TABLE IF EXISTS $production_name"."_32");
$search_dbh->do("DROP TABLE IF EXISTS $production_name"."_255");
$search_dbh->do("DELETE FROM multi WHERE production_name = ".$search_dbh->quote($production_name));
$search_dbh->do("DELETE FROM multi_32 WHERE production_name = ".$search_dbh->quote($production_name));
$search_dbh->do("DELETE FROM multi_255 WHERE production_name = ".$search_dbh->quote($production_name));


# loop through seq regions and genes adding search data to search db
$search_dbh->do("CREATE TABLE $production_name LIKE multi");
$search_dbh->do("CREATE TABLE $production_name"."_32 LIKE multi_32");
$search_dbh->do("CREATE TABLE $production_name"."_255 LIKE multi_255");
index_meta($core_dbh,$search_dbh,$production_name);
index_seq_regions($core_dbh,$search_dbh,$production_name);
index_genes($core_dbh,$search_dbh,$production_name);
