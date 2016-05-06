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
  'FILES' => 	{	'GFF' => 1,
          'SCAFFOLD' => [ 'CONTIG' ],
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

# split gff to extract FASTA sequence if necessary
if ($params{'GFF'}{'SPLIT'}){
	split_gff($params,\%infiles);
}


## generate summary of gff file
foreach my $file (keys %infiles){
	if ($infiles{$file}{'type'} eq 'gff'){
		## generate feature summary and write to file
		print STDERR "Calculating summary statistics on [FILES] $file $infiles{$file}{'name'}\n";
		gff_feature_summary($params,$infiles{'GFF'}->{'name'});
	}
}



sub usage {
	return "USAGE: perl /path/to/summarise_files.pl /pat/to/config_file.ini";
}
