#!/usr/bin/perl -w

use strict;
use Cwd 'abs_path';
use File::Basename;
use JSON;

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

  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
my %infiles;
my $features = undef;
my @inis = @ARGV;
my %stats;
mkdir "web";
@ARGV = ();
while (my $ini_file = shift @inis){
	load_ini($params,$ini_file,\%sections,scalar(@ARGV));
	## download/obtain files using methods suggested by file paths and extensions
  foreach my $subsection (sort keys %{$params->{'FILES'}}){
    if ($subsection eq 'GFF' || $subsection eq 'CEGMA' || $subsection eq 'BUSCO'){
 	    ($infiles{$subsection}{'name'},$infiles{$subsection}{'type'}) = fetch_file($params->{'FILES'}{$subsection});
    }
  }
}

my $production_name = $params->{'META'}{'SPECIES.PRODUCTION_NAME'};
my $exportdir = 'exported';


my $dbh = core_db_connect($params);
my $sth = $dbh->prepare('SELECT count(*) FROM gene');
$sth->execute;
$stats{'gene_count'} = $sth->fetchrow_arrayref->[0];
$sth = $dbh->prepare('SELECT count(*) FROM transcript');
$sth->execute;
$stats{'transcript_count'} = $sth->fetchrow_arrayref->[0];
$sth = $dbh->prepare('SELECT count(*) FROM translation');
$sth->execute;
$stats{'translation_count'} = $sth->fetchrow_arrayref->[0];

## generate summaries of CEGMA files
my $cegma;
foreach my $file (keys %infiles){
	# 1.1: Read CEGMA scores from completness_report
	if ($file =~ m/^CEGMA$/ && $infiles{$file}{'type'} eq 'txt'){
		$cegma = cegma_file_summary($infiles{$file}->{'name'});
		foreach my $key (keys %{$cegma}){
			$stats{$key} = $cegma->{$key};
		}
	}
}


$infiles{'SCAFFOLD'}->{'name'} = "$exportdir/$production_name"."_-_scaffolds.fa";

if (-s $infiles{'SCAFFOLD'}->{'name'}){
	## calculate scaffold stats and write to file
	print STDERR "Calculating summary statistics on $infiles{'SCAFFOLD'}{'name'}\n";
	my $tmp_stats = fasta_file_summary($params,$infiles{'SCAFFOLD'},'SCAFFOLD',$cegma);
  my $json = JSON->new;
  $json->pretty(1);

  open JS,">web/$production_name.stats.json";
  print JS $json->encode($tmp_stats),"\n";
  close JS;
	foreach my $key (keys %{$tmp_stats}){
		$stats{$key} = $tmp_stats->{$key};
	}
}


## generate template about_production_name.html page
about_page($params);

## generate stats_production_name.html page - IPTop500 may not exist yet but create link anyway
stats_page($params,\%stats);

## generate production_name_assembly.html page with assembly visualisation (and description)
# 1.1: Added GC wiggle to assembly badge
# 1.1: Added gff feature histogram to html
assembly_page($params,\%stats);

# 1.1: Generate ini file for lepbase-ensembl
web_ini_file($params,\%stats);

sub usage {
	return "USAGE: perl /path/to/generate_file_stats.pl /pat/to/config_file.ini";
}
