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
        load_ini($params,$ini_file,\%sections,scalar(@ARGV));
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
    -user   => $params->{'DATABASE_CORE'}{'RW_USER'},
    -pass   => $params->{'DATABASE_CORE'}{'RW_PASS'},
    -dbname => $dbname,
    -host   => $params->{'DATABASE_CORE'}{'HOST'},
    -port   => $params->{'DATABASE_CORE'}{'PORT'},
    -driver => 'mysql'
);


my $meta_container = $dba->get_adaptor("MetaContainer");

my $display_name    = $meta_container->get_display_name();
my $assembly_name   = $meta_container->single_value_by_key('ASSEMBLY.NAME');
$display_name .= ' '.$assembly_name;

# convert display name spaces to underscores
# $display_name =~ s/ /_/g;

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
      push @{$features{'Transcripts'}->{'lengths'}},$transcript->length();
      $features{'Transcripts'}->{'base_count'} = base_composition($transcript->seq->seq(),$features{'Transcripts'}->{'base_count'});
      my $translateable_seq;
      if ($translateable_seq = $transcript->translateable_seq()){
        $features{'codon_count'} = codon_count($translateable_seq,$features{'codon_count'});
        push @{$features{'CDS'}->{'lengths'}},length($translateable_seq);
        $features{'CDS'}->{'base_count'} = base_composition($translateable_seq,$features{'CDS'}->{'base_count'});
      }
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

my %cegma_busco;
if ($meta_container->single_value_by_key('assembly.busco_complete')){
  my $busco;
  $busco->{'C'} = $meta_container->single_value_by_key('assembly.busco_complete')*1;
  $busco->{'D'} = $meta_container->single_value_by_key('assembly.busco_duplicated')*1;
  $busco->{'F'} = $meta_container->single_value_by_key('assembly.busco_fragmented')*1;
  $busco->{'M'} = $meta_container->single_value_by_key('assembly.busco_missing')*1;
  $busco->{'n'} = $meta_container->single_value_by_key('assembly.busco_number')*1;
  $cegma_busco{'busco'} = $busco;
}
elsif ($meta_container->single_value_by_key('assembly.cegma_complete')){
  $cegma_busco{'cegma_complete'} = $meta_container->single_value_by_key('assembly.cegma_complete')*1;
  $cegma_busco{'cegma_partial'} = $meta_container->single_value_by_key('assembly.cegma_partial')*1;
}
my $assembly_stats = scaffold_summary($params,\@scaffolds,'SCAFFOLD',\%cegma_busco);
$assembly_stats->{'name'} = $display_name;

my $json = JSON->new;
$json->pretty(1);

mkdir 'web';
open JS,">web/$dbname.assembly-stats.json";
print JS $json->encode($assembly_stats),"\n";
close JS;

# add some assembly stats to meta table
$meta_container->store_key_value('assembly.span',$assembly_stats->{'assembly'});
$meta_container->store_key_value('assembly.gc_percent',$assembly_stats->{'GC'});
$meta_container->store_key_value('assembly.atgc',$assembly_stats->{'ATGC'});
$meta_container->store_key_value('assembly.n',$assembly_stats->{'N'});
$meta_container->store_key_value('assembly.scaffold_count',scalar @scaffolds);

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
  for my $l (@{$features{$key}->{'lengths'}}) {
    my $bin = ceil($binres*log10($l > 0 ? $l : 1))/$binres;
    $bin = FormatSigFigs(10**$bin,1);
    $bin =~ s/\.$//;
    $features{$key}->{'binned'}[$bins{$bin}]++;
  }
  delete $features{$key}->{'lengths'};

  print $maxbin,"\n";
  print 10**$maxbin,"\n";
}
$features{'name'} = $display_name;
$json = JSON->new;
$json->pretty(1);

open JS,">web/$dbname.codon-usage.json";
print JS $json->encode(\%features),"\n";
close JS;


$meta_container->store_key_value('genebuild.gene_count',sum @{$features{'Genes'}->{'binned'}});
$meta_container->store_key_value('genebuild.transcript_count',sum @{$features{'Transcripts'}->{'binned'}});
$meta_container->store_key_value('genebuild.cds_count',sum @{$features{'CDS'}->{'binned'}});
$meta_container->store_key_value('genebuild.exon_count',sum @{$features{'Exons'}->{'binned'}});

my %meta;

$meta{'provider'}->{'name'} = $meta_container->single_value_by_key('provider.name');
$meta{'provider'}->{'url'} = $meta_container->single_value_by_key('provider.url');
$meta{'species'}->{'common_name'} = $meta_container->single_value_by_key('species.common_name');
$meta{'species'}->{'display_name'} = $meta_container->single_value_by_key('species.display_name');
$meta{'species'}->{'scientific_name'} = $meta_container->single_value_by_key('species.scientific_name');
$meta{'species'}->{'taxonomy_id'} = $meta_container->single_value_by_key('species.taxonomy_id');
$meta{'assembly'}->{'accession'} = $meta_container->single_value_by_key('assembly.accession');
$meta{'assembly'}->{'date'} = $meta_container->single_value_by_key('assembly.date');
$meta{'assembly'}->{'name'} = $meta_container->single_value_by_key('assembly.name');
$meta{'assembly'}->{'span'} = $meta_container->single_value_by_key('assembly.span');
$meta{'assembly'}->{'gc_percent'} = $meta_container->single_value_by_key('assembly.gc_percent');
$meta{'assembly'}->{'n'} = $meta_container->single_value_by_key('assembly.n');
$meta{'assembly'}->{'atgc'} = $meta_container->single_value_by_key('assembly.atgc');
$meta{'assembly'}->{'scaffold_count'} = $meta_container->single_value_by_key('assembly.scaffold_count');
if ($meta_container->single_value_by_key('assembly.cegma_complete')){
  $meta{'assembly'}->{'cegma_complete'} = $meta_container->single_value_by_key('assembly.cegma_complete');
  $meta{'assembly'}->{'cegma_partial'} = $meta_container->single_value_by_key('assembly.cegma_partial');
}
if ($meta_container->single_value_by_key('assembly.busco_complete')){
  $meta{'assembly'}->{'busco_complete'} = $meta_container->single_value_by_key('assembly.busco_complete');
  $meta{'assembly'}->{'busco_duplicated'} = $meta_container->single_value_by_key('assembly.busco_duplicated');
  $meta{'assembly'}->{'busco_fragmented'} = $meta_container->single_value_by_key('assembly.busco_fragmented');
  $meta{'assembly'}->{'busco_missing'} = $meta_container->single_value_by_key('assembly.busco_missing');
  $meta{'assembly'}->{'busco_number'} = $meta_container->single_value_by_key('assembly.busco_number');
}
$meta{'genebuild'}->{'method'} = $meta_container->single_value_by_key('genebuild.method');
$meta{'genebuild'}->{'start_date'} = $meta_container->single_value_by_key('genebuild.start_date');
$meta{'genebuild'}->{'version'} = $meta_container->single_value_by_key('genebuild.version');
$meta{'genebuild'}->{'gene_count'} = $meta_container->single_value_by_key('genebuild.gene_count');
$meta{'genebuild'}->{'transcript_count'} = $meta_container->single_value_by_key('genebuild.transcript_count');
$meta{'genebuild'}->{'cds_count'} = $meta_container->single_value_by_key('genebuild.cds_count');
$meta{'genebuild'}->{'exon_count'} = $meta_container->single_value_by_key('genebuild.exon_count');

$json = JSON->new;
$json->pretty(1);

open JS,">web/$dbname.meta.json";
print JS $json->encode(\%meta),"\n";
close JS;

# interpro top 500 hits
my $dbh = core_db_connect($params);
my $sth = $dbh->prepare('SELECT COUNT(*) AS count, ipr.name, ipr.description
                         FROM (
                           SELECT DISTINCT t.gene_id AS gene_id, x.dbprimary_acc AS name, x.description AS description
                           FROM transcript AS t
                           JOIN object_xref AS ox ON ox.ensembl_id = t.transcript_id
                           JOIN xref AS x ON ox.xref_id = x.xref_id
                           WHERE x.dbprimary_acc LIKE \'IPR%\'
                         ) AS ipr
                         GROUP BY ipr.name
                         ORDER BY count DESC
                         LIMIT 500');
$sth->execute();
if ($sth->rows > 0){
  my @ipr;
  while (my $row = $sth->fetchrow_arrayref()){
    push @ipr,{count => $row->[0]*1,name => $row->[1],description => $row->[2]};
  }
  $json = JSON->new;
  $json->pretty(1);

  open JS,">web/$display_name.IPtop500.json";
  print JS $json->encode(\@ipr),"\n";
  close JS;
}

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
