#!/usr/bin/perl -w

use strict;
use Cwd 'abs_path';
use File::Basename;
use Module::Load;
use Parallel::ForkManager;

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
  'DATABASE_CORE' =>	{ 	'HOST' => 1,
              'PORT' => 1,
              'RO_USER' => 1
            },
   'TAXA' =>	{},
   'SETUP' =>	{  'FASTA_DIR' => 1,
                 'REMOVE' =>	1
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
push @INC, $lib;
load Bio::EnsEMBL::DBSQL::DBAdaptor;


my $outdir = $params->{'SETUP'}{'FASTA_DIR'};
mkdir $outdir unless -d $outdir;

my %overwrite;
for my $taxon (@{$params->{'SETUP'}{'REMOVE'}}){
  $overwrite{$taxon} = 1;
}

# Max 2 processes for parallel download
my $pm = new Parallel::ForkManager(2);

for my $taxon (keys %{$params->{'TAXA'}}){

  $pm->start and next; # do the fork

  # test if file exists and should be overwritten
  if (-e "$outdir/$taxon\_-_canonical_proteins.fa"){
    next unless $overwrite{$taxon};
  }

  my $dbname = $params->{'TAXA'}{$taxon};

  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
      -user   => $params->{'DATABASE_CORE'}{'RO_USER'},
      -dbname => $dbname,
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

  # Get all transcripts
  # For each transcript
  #   Print to proteins.fa, cds.fa, cds_translationid.fa

  my $transcript_adaptor = $dba->get_TranscriptAdaptor();
  my @transcripts        = @{$transcript_adaptor->fetch_all_by_biotype('protein_coding')};
  my ($pep, $cds, $bounded_exon, $transcript_id, $translation_id, $desc)  = ("","","","","");
  my ($canonical_protein_fh, $bounded_exon_fh, $cds_translationid_fh);
  my $canonical_count = 0;

  system "mkdir -p $outdir/canonical_proteins $outdir/canonical_cds_translationid $outdir/canonical_protein_bounded_exon";
  open $canonical_protein_fh, ">", "$outdir/canonical_proteins/$taxon\_-_canonical_proteins.fa"   or die $!;
  open $cds_translationid_fh, ">", "$outdir/canonical_cds_translationid/$taxon\_-_canonical_cds_translationid.fa"    or die $!;
  open $bounded_exon_fh,      ">", "$outdir/canonical_protein_bounded_exon/$taxon\_-_canonical_protein_bounded_exon.fa" or die $!;

  foreach my $transcript (@transcripts) {
    if (defined $transcript->translate() ) {
      if ($transcript->is_canonical()){
        $transcript_id   = $transcript->stable_id();
        $translation_id  = $transcript->translation()->stable_id();
        $pep = $transcript->translate()->seq;
        $cds = $transcript->translateable_seq();
        $bounded_exon = _prepare_exon_sequences($transcript,$pep);
        print $canonical_protein_fh ">${taxon}_$translation_id $dbname protein\n$pep\n";
        print $cds_translationid_fh ">${taxon}_$translation_id $dbname cds_translationid\n$cds\n";
        print $bounded_exon_fh      ">${taxon}_$translation_id $dbname protein_bounded_exon\n$bounded_exon\n";
        $canonical_count++;
      }
    }
  }
  print "$dbname - Num of canonical proteins : $canonical_count\n";

  $pm->finish; # do the exit in the child process
}

$pm->wait_all_children;

sub usage {
	return "USAGE: perl /path/to/export_core_sequences.pl /path/to/config_file.ini";
}

sub _prepare_exon_sequences {


#    my $self = shift;
#
#    # If there is the exon_bounded sequence, it is only a matter of splitting it and alternating the case
#    my $exon_bounded_seq = $self->{_sequence_exon_bounded};
#    $exon_bounded_seq = $self->adaptor->db->get_SequenceAdaptor->fetch_other_sequence_by_member_id_type($self->seq_member_id, 'exon_bounded') unless $exon_bounded_seq;
#
#    if ($exon_bounded_seq) {
#        $self->{_sequence_exon_bounded} = $exon_bounded_seq;
#        my $i = 0;
#        $self->{_sequence_exon_cased} = join('', map {$i++%2 ? lc($_) : $_} split( /[boj]/, $exon_bounded_seq));
#
#    } else {
#
#        my $sequence = $self->sequence;
        my ($transcript,$sequence) = @_;
        my $sequence_exon_bounded;
        my @exons = @{$transcript->get_all_translateable_Exons};
        # @exons probably doesn't match the protein if there are such edits
        my @seq_edits = @{$transcript->translation->get_all_SeqEdits('amino_acid_sub')};
        push @seq_edits, @{$transcript->get_all_SeqEdits('_rna_edit')};

        if (((scalar @exons) <= 1) or (scalar(@seq_edits) > 0)) {
            $sequence_exon_bounded = $sequence;
            return $sequence_exon_bounded;
        }

        # Otherwise, we have to parse the exons
        my %boundary_chars = (0 => 'o', 1 => 'b', 2 => 'j');
        my $left_over = $exons[0]->phase > 0 ? -$exons[0]->phase : 0;
        my @this_seq = ();
       # my @exon_sequences = ();
        foreach my $exon (@exons) {
            my $exon_pep_len = POSIX::ceil(($exon->length - $left_over) / 3);
            my $exon_seq = substr($sequence, 0, $exon_pep_len, '');
            $left_over += 3*$exon_pep_len - $exon->length;
            #printf("%s: exon of len %d -> phase %d: %s\n", $transcript->stable_id, $exon_pep_len, $left_over, $exon_seq);
            push @this_seq, $exon_seq;
            push @this_seq, $boundary_chars{$left_over};
           # push @exon_sequences, scalar(@exon_sequences)%2 ? $exon_seq : lc($exon_seq);
            die sprintf('Invalid phase: %s', $left_over) unless exists $boundary_chars{$left_over}
        }
        die sprintf('%d characters left in the sequence of %s', length($sequence), $transcript->stable_id) if $sequence;
        pop @this_seq;
        $sequence_exon_bounded = join('', @this_seq);
        #$sequence_exon_cased = join('', @exon_sequences);

}
