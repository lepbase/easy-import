#!/usr/bin/perl -w

use strict;
use Cwd 'abs_path';
use File::Basename;
use List::Util qw(max);
use POSIX qw(ceil);
use JSON;
use Math::SigFigs;
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
  'ENSEMBL' =>  {       'LOCAL' => 1
          },
  'DATABASE_CORE' =>    {       'NAME' => 1,
              'HOST' => 1,
              'PORT' => 1,
              'RW_USER' => 1,
              'RW_PASS' => 1,
              'RO_USER' => 1
            }
  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
        load_ini($params,$ini_file,\%sections);
}


my $lib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl/modules';
my $iolib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-io/modules';
#my $comparalib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-compara/modules';
push @INC, $lib;
push @INC, $iolib;
#push @INC, $comparalib;
load Bio::EnsEMBL::Registry;
load Bio::EnsEMBL::DBSQL::DBAdaptor;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => $params->{'DATABASE_CORE'}{'HOST'},
    -user => $params->{'DATABASE_CORE'}{'RO_USER'},
    -port   => $params->{'DATABASE_CORE'}{'PORT'},
    -driver => 'mysql',
);

my $dbname = $params->{'DATABASE_CORE'}{'NAME'};

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => $params->{'DATABASE_CORE'}{'RO_USER'},
    -dbname => $dbname,
    -host   => $params->{'DATABASE_CORE'}{'HOST'},
    -port   => $params->{'DATABASE_CORE'}{'PORT'},
    -driver => 'mysql'
);


my $meta_container = $dba->get_adaptor("MetaContainer");

my $display_name    = $meta_container->get_display_name();

# convert display name spaces to underscores
$display_name =~ s/ /_/g;

my $ontology_adaptor = $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
my $slice_adaptor    = $dba->get_adaptor("Slice");
my $supercontigs     = $slice_adaptor->fetch_all('toplevel');
my @scaffolds;
my %features;
$features{'codon_count'} = {};


foreach my $slice (@{$supercontigs}) {
  push @scaffolds,$slice->seq();
  push @{$features{'Scaffolds'}->{'lengths'}},$slice->end();
  $features{'Scaffolds'}->{'base_count'} = base_composition($slice->seq(),$features{'Scaffolds'}->{'base_count'});
  my $genes = $slice->get_all_Genes;
  while ( my $gene = shift @{$genes} ) {
    push @{$features{'Genes'}->{'lengths'}},$gene->length();
    $features{'Genes'}->{'base_count'} = base_composition($gene->seq(),$features{'Genes'}->{'base_count'});
    my $transcripts = $gene->get_all_Transcripts();
    while ( my $transcript = shift @{$transcripts} ) {
      my $translateable_seq = $transcript->translateable_seq();
      push @{$features{'Transcripts'}->{'lengths'}},$transcript->length();
      $features{'Transcripts'}->{'base_count'} = base_composition($transcript->seq->seq(),$features{'Transcripts'}->{'base_count'});
      $features{'codon_count'} = codon_count($translateable_seq,$features{'codon_count'});
      push @{$features{'CDS'}->{'lengths'}},length($translateable_seq);
      $features{'CDS'}->{'base_count'} = base_composition($translateable_seq,$features{'CDS'}->{'base_count'});
      foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
        push @{$features{'Exons'}->{'lengths'}},$exon->length();
        $features{'Exons'}->{'base_count'} = base_composition($exon->seq->seq(),$features{'Exons'}->{'base_count'});
      }
    }
  }
}
$features{'Scaffolds'}->{'base_count'}->{'A'} = $features{'Scaffolds'}->{'base_count'}->{'A'} + $features{'Scaffolds'}->{'base_count'}->{'U'};
$features{'Scaffolds'}->{'base_count'}->{'C'} = $features{'Scaffolds'}->{'base_count'}->{'C'} + $features{'Scaffolds'}->{'base_count'}->{'G'};
$features{'Scaffolds'}->{'base_count'}->{'G'} = $features{'Scaffolds'}->{'base_count'}->{'C'};
$features{'Scaffolds'}->{'base_count'}->{'U'} = $features{'Scaffolds'}->{'base_count'}->{'A'};

my $assembly_stats = scaffold_summary($params,\@scaffolds,'SCAFFOLD');#,$cegma);
my $json = JSON->new;
$json->pretty(1);

mkdir 'web';
open JS,">web/$display_name.assembly-stats.json";
print JS $json->encode($assembly_stats),"\n";
close JS;

foreach my $key (keys %features){
  next if $key eq 'codon_count';
  #my $sdls = Statistics::Descriptive::LogScale->new ();
  #$sdls->add_data(@{$features{$key}{'lengths'}});
  my $max = max(@{$features{$key}->{'lengths'}});
  print $key,"\n";
  print $max,"\n";
  my $maxbin = ceil(log10($max));
  my $binres = $maxbin <= 7 ? $maxbin <= 4 ? 3 : 2 : 1;
  my $binsize = 1 / $binres;
  my %bins;
  my $i = 0;
  for (my $b = $binsize; $b <= $maxbin; $b += $binsize){
    my $bin = FormatSigFigs(10**$b,1);
    $bin =~ s/\.$//;
    push @{$features{$key}->{'bins'}},$bin;
    $features{$key}->{'binned'}[$i] = 0;
    $bins{$bin} = $i;
    $i++;
  }
  while (my $l = shift @{$features{$key}->{'lengths'}}){
    my $bin = ceil($binres*log10($l))/$binres;
    $bin = FormatSigFigs(10**$bin,1);
    $bin =~ s/\.$//;
    $features{$key}->{'binned'}[$bins{$bin}]++;
  }
  delete $features{$key}->{'lengths'};

  print $maxbin,"\n";
  print 10**$maxbin,"\n";
}

$json = JSON->new;
$json->pretty(1);

open JS,">web/$display_name.codon-usage.json";
print JS $json->encode(\%features),"\n";
close JS;

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

sub base_composition {
  my $str = shift;
  my $bases = shift || {};
  $bases->{'A'} += () = $str =~ /a/gi;
  $bases->{'C'} += () = $str =~ /c/gi;
  $bases->{'G'} += () = $str =~ /g/gi;
  $bases->{'U'} += () = $str =~ /[ut]/gi;
  return $bases;
}

sub codon_count {
  my $str = shift;
  my $codons = shift;
  $str =~ s/T/U/ig;
  while ($str =~ s/^([ACGU]{3})//i){
    $codons->{$1}++;
  }
  return $codons;
}

__END__
mkdir 'exported';
my   $output_fh;
open $output_fh, ">exported/$display_name.gff";

my   $ontology_adaptor = $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
my   $serializer       = Bio::EnsEMBL::Utils::IO::GFFSerializer->new($ontology_adaptor,$output_fh);

my   $slice_adaptor    = $dba->get_adaptor("Slice");
my   $supercontigs     = $slice_adaptor->fetch_all('toplevel');

foreach my $slice (@{$supercontigs}) {
    my $genes = $slice->get_all_Genes;
    while ( my $gene = shift @{$genes} ) {
        $serializer->print_feature($gene);
        my $transcripts = $gene->get_all_Transcripts();
        while ( my $transcript = shift @{$transcripts} ) {
            $serializer->print_feature($transcript);
            foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
                $serializer->print_feature($exon);
            }
            foreach my $cds ( @{ $transcript->get_all_CDS() } ) {
                 $serializer->print_feature($cds);
            }
        }
    }
}


sub usage {
	return "USAGE: perl /path/to/export_gff.pl /path/to/config_file.ini";
}
