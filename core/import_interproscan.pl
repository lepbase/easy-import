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
use EasyImport::Xref;

## load parameters from an INI-style config file
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file);
}

## download/obtain files using methods suggested by file paths and extensions
my %infiles;
foreach my $subsection (sort keys %{$params->{'FILES'}}){
	($infiles{$subsection}{'name'},$infiles{$subsection}{'type'}) = fetch_file($params->{'FILES'}{$subsection});
}

## check that an IPRSCAN file has been specified
die "ERROR: an interproscan result file must be specified at [FILES] IPRSCAN\n" unless $infiles{'IPRSCAN'};

## connect to the db
my $dbh = core_db_connect($params);

## load the IPRSCAN file into the database
push @ARGV,$infiles{'IPRSCAN'}{'name'};
my $hits = read_iprscan($dbh,$params);

if ($hits){
	## print IPR top 500 table
	make_IPtop500_table($hits,$params);
}





__DATA__
nBa.0.1-t07619-RA	3cf6e457984e03aca1005f66feb42f94	313	Pfam	PF00151	Lipase	36	274	1.2E-35	T	25-02-2015
nBa.0.1-t07619-RA	3cf6e457984e03aca1005f66feb42f94	313	PRINTS	PR00821	Triacylglycerol lipase family signature	41	60	8.8E-12	T	25-02-2015
nBa.0.1-t07619-RA	3cf6e457984e03aca1005f66feb42f94	313	PRINTS	PR00821	Triacylglycerol lipase family signature	242	257	8.8E-12	T	25-02-2015
nBa.0.1-t07619-RA	3cf6e457984e03aca1005f66feb42f94	313	PRINTS	PR00821	Triacylglycerol lipase family signature	135	153	8.8E-12	T	25-02-2015
nBa.0.1-t07619-RA	3cf6e457984e03aca1005f66feb42f94	313	SUPERFAMILY	SSF53474		22	274	2.76E-34	T	25-02-2015
nBa.0.1-t07619-RA	3cf6e457984e03aca1005f66feb42f94	313	Gene3D	G3DSA:3.40.50.1820		21	307	1.1E-55	T	25-02-2015
nBa.0.1-t02120-RC	ed294fc1dc72020193a7c1758bca6369	889	TMHMM	TMhelix	Region of a membrane-bound protein predicted to be embedded in the membrane.	584	602	-	T	25-02-2015
nBa.0.1-t02120-RC	ed294fc1dc72020193a7c1758bca6369	889	TMHMM	TMhelix	Region of a membrane-bound protein predicted to be embedded in the membrane.	544	563	-	T	25-02-2015
nBa.0.1-t02120-RC	ed294fc1dc72020193a7c1758bca6369	889	Phobius	TRANSMEMBRANE	Region of a membrane-bound protein predicted to be embedded in the membrane.	746	766	-	T	25-02-2015
nBa.0.1-t02120-RC	ed294fc1dc72020193a7c1758bca6369	889	ProSiteProfiles	PS50893	ATP-binding cassette, ABC transporter-type domain profile.	76	323	16.367	T	25-02-2015
