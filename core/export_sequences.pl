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
  'META' =>	{ 	'SPECIES.PRODUCTION_NAME' => 1
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
foreach my $subsection (sort keys %{$params->{'FILES'}}){
	($infiles{$subsection}{'name'},$infiles{$subsection}{'type'}) = fetch_file($params->{'FILES'}{$subsection});
}

use lib $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl/modules';
use Bio::EnsEMBL::DBSQL::DBAdaptor;

=head1

extract_blast_sequences.pl

Given an ensembl species database name, extracts 4 files named after Species_display_name
with headers formatted as:

Species_display_name_-_scaffolds.fa
>scaffoldname dbname scaffold

Species_display_name_-_proteins.fa
>translationid dbname protein description with spaces

Species_display_name_-_cds.fa
>transcriptid dbname cds description with spaces

Species_display_name_-_cds_translationid.fa
>translationid dbname cds_translationid description with spaces

Notes:

1. Advantage of using dbname is that we get the production name AND genebuild in the header
2. Need for cds_translationid.fa is that some orthology pipelines expects the SAME id for protein and cds

USAGE:

extract_blast_sequences.pl heliconius_melpomene_hmel2_27_80_1

# note for ruby script - the seq region in the blast results must have start < end for external URL link to work

=cut

my $dbname = shift @ARGV;

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => $params->{'DATABASE_CORE'}{'RO_USER'},
    -dbname => lc $params->{'META'}{'SPECIES.PRODUCTION_NAME'},
    -host   => $params->{'DATABASE_CORE'}{'HOST'},
    -port   => $params->{'DATABASE_CORE'}{'PORT'},
    -driver => 'mysql'
);

my $meta_container = $dba->get_adaptor("MetaContainer");

# Get meta names from db

my $production_name = $meta_container->get_production_name();
my $scientific_name = $meta_container->get_scientific_name();
my $display_name    = $meta_container->get_display_name();

# convert display name spaces to underscores
$display_name =~ s/ /_/g;

# Get all scaffolds

my $slice_adaptor = $dba->get_SliceAdaptor();
my @supercontigs  = @{$slice_adaptor->fetch_all('toplevel')};
my $supercontig_count = 0;

open (SCAFFOLDS, ">", "$display_name\_-_scaffolds.fa") or die $!;

foreach my $slice (@supercontigs) {
    print SCAFFOLDS ">" . $slice->seq_region_name() . " $dbname scaffold\n" . $slice->seq() . "\n";
    $supercontig_count++;
}

print "$dbname - Num of supercontigs : $supercontig_count\n";

# Get all transcripts
# For each transcript
#   Ignore source (we may want to come back to this if there is a need for separate blast databases by source)
#   Print to proteins.fa, cds.fa, cds_translationid.fa

my $transcript_adaptor = $dba->get_TranscriptAdaptor();
my $gene_adaptor       = $dba->get_GeneAdaptor();
my @transcripts        = @{$transcript_adaptor->fetch_all_by_biotype('protein_coding')};
my $gene;
my ($pep, $cds, $transcript_id, $translation_id, $desc)  = ("","","","","");
my ($protein_fh, $cds_fh, $cds_translationid_fh);
my $protein_count = 0;

open $protein_fh,           ">", "$display_name\_-_proteins.fa"          or die $!;
open $cds_fh,               ">", "$display_name\_-_cds.fa"               or die $!;
open $cds_translationid_fh, ">", "$display_name\_-_cds_translationid.fa" or die $!;

foreach my $transcript (@transcripts) {
    if (defined $transcript->translate() ) {
        $transcript_id   = $transcript->stable_id();
        $translation_id  = $transcript->translation()->stable_id();
        $desc = "";
        # Leave desc out for now, getting a stack overflow error on next line
        # $gene            = $gene_adaptor->fetch_by_transcript_id($transcript_id);
        # if    (defined $transcript->description) {
        #        $desc = $transcript->description;
        # }
        # elsif (defined $gene->description) {
        #        $desc = $gene->description;
        # }
        $pep = $transcript->translate()->seq;
        $cds = $transcript->translateable_seq();
        # print $cds_fh               ">$transcript_id $dbname cds $desc\n$cds\n";
        # print $cds_translationid_fh ">$translation_id $dbname cds_translationid $desc\n$cds\n";
        # print $protein_fh           ">$translation_id $dbname protein $desc\n$pep\n";
        print $cds_fh               ">$transcript_id $dbname cds $desc\n$cds\n";
        print $cds_translationid_fh ">$translation_id $dbname cds_translationid $desc\n$cds\n";
        print $protein_fh           ">$translation_id $dbname protein $desc\n$pep\n";
        $protein_count++;
    }
}

print "$dbname - Num of proteins     : $protein_count\n";


sub usage {
	return "USAGE: perl /path/to/export_sequences.pl /path/to/config_file.ini";
}
