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

my $lib = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl/modules';
push @INC, $lib;
load Bio::EnsEMBL::DBSQL::DBAdaptor;

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

my $dbname = $params->{'DATABASE_CORE'}{'NAME'};

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

my $outdir = 'exported';
mkdir $outdir;

# Get all scaffolds

my $slice_adaptor = $dba->get_SliceAdaptor();
my @supercontigs  = @{$slice_adaptor->fetch_all('toplevel')};
my $supercontig_count = 0;

open (SCAFFOLDS, ">", "$outdir/$display_name\_-_scaffolds.fa") or die $!;

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
my ($pep, $cds, $bounded_exon, $transcript_id, $translation_id, $desc)  = ("","","","","");
my ($protein_fh, $cds_fh, $bounded_exon_fh, $cds_translationid_fh);
my $protein_count = 0;

open $protein_fh,           ">", "$outdir/$display_name\_-_proteins.fa"          or die $!;
open $cds_fh,               ">", "$outdir/$display_name\_-_cds.fa"               or die $!;
open $cds_translationid_fh, ">", "$outdir/$display_name\_-_cds_translationid.fa" or die $!;
open $bounded_exon_fh,      ">", "$outdir/$display_name\_-_protein_bounded_exon.fa"  or die $!;

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
        $bounded_exon = _prepare_exon_sequences($transcript,$pep);
        # print $cds_fh               ">$transcript_id $dbname cds $desc\n$cds\n";
        # print $cds_translationid_fh ">$translation_id $dbname cds_translationid $desc\n$cds\n";
        # print $protein_fh           ">$translation_id $dbname protein $desc\n$pep\n";
        print $cds_fh               ">$transcript_id $dbname cds $desc\n$cds\n";
        print $cds_translationid_fh ">$translation_id $dbname cds_translationid $desc\n$cds\n";
        print $bounded_exon_fh      ">$translation_id $dbname protein_bounded_exon $desc\n$bounded_exon\n";
        print $protein_fh           ">$translation_id $dbname protein $desc\n$pep\n";
        $protein_count++;
    }
}

print "$dbname - Num of proteins     : $protein_count\n";


sub usage {
	return "USAGE: perl /path/to/export_sequences.pl /path/to/config_file.ini";
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
