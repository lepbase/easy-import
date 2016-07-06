#!/usr/bin/perl -w

use strict;
use Cwd 'abs_path';
use File::Basename;
use Module::Load;
use Bio::SeqIO;
use Text::LevenshteinXS qw(distance);


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
  'FILES' =>	{ 	'PROTEIN' => 1
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

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => $params->{'DATABASE_CORE'}{'RW_USER'},
    -pass   => $params->{'DATABASE_CORE'}{'RW_PASS'},
    -dbname => $params->{'DATABASE_CORE'}{'NAME'},
    -host   => $params->{'DATABASE_CORE'}{'HOST'},
    -port   => $params->{'DATABASE_CORE'}{'PORT'},
    -driver => 'mysql'
);

my $meta_container = $dba->get_adaptor("MetaContainer");

my $display_name    = $meta_container->get_display_name();
$display_name =~ s/ /_/g;

## download/obtain files using methods suggested by file paths and extensions
my %infiles;
($infiles{'PROTEIN'}{'name'},$infiles{'PROTEIN'}{'type'}) = fetch_file($params->{'FILES'}{'PROTEIN'});

my $provider = $infiles{'PROTEIN'}{'name'};

####

#print LOG "\n----\nChecking for duplicates in $provider\n----\n";

my $providerin = Bio::SeqIO->new(
                            -file   => "<$provider",
                            -format => "FASTA",
                            );

my %providerhash;
my @providerdups;

while ( my $seq = $providerin->next_seq) {
    my $id  = $seq->display_id;
    my $seqstr = $seq->seq;
    $seqstr =~ s/\W*$//;
    if (exists $providerhash{$id}) {
        print LOG "DUPLICATE $id\n";
        push  @providerdups, $id;
    }
    else {
        $providerhash{$id} = $seqstr;
    }
}

my $outdir = 'summary';
mkdir $outdir;

open  LOG, ">$outdir/$display_name"."_-_proteins.log" or die $!;

print LOG "Num of unique     seq IDs in $provider: " . scalar (keys %providerhash) . "\n";
print LOG "Num of duplicated seq IDs in $provider: " . scalar (@providerdups) . "\n";

if (scalar @providerdups > 0) {
    open  DUP1, ">$outdir/$provider.dups.ids" or die $!;
    print DUP1  join("\n",@providerdups)  . "\n";
}

####

#print LOG "\n----\nChecking for duplicates in $exported\n----\n";

my %exportedhash;
my @exporteddups;

my $transcript_adaptor = $dba->get_TranscriptAdaptor();
my @transcripts        = @{$transcript_adaptor->fetch_all_by_biotype('protein_coding')};
foreach my $transcript (@transcripts) {
  if (defined $transcript->translate() ) {
    my $id = $transcript->translation()->stable_id();
    my $seqstr = $transcript->translate()->seq;
    $seqstr =~ s/\W*$//;
    if (exists $exportedhash{$id}) {
      print LOG "DUPLICATE $id\n";
      push  @exporteddups, $id;
    }
    else {
      $exportedhash{$id} = $seqstr;
    }
  }
}

print LOG "Num of unique     seq IDs in database: " . scalar (keys %exportedhash) . "\n";
print LOG "Num of duplicated seq IDs in database: " . scalar (@exporteddups) . "\n";

if (scalar @exporteddups > 0) {
    open  DUP2, ">$outdir/database.dups.ids" or die $!;
    print DUP2  join("\n",@exporteddups)  . "\n";
}

####

#print LOG "\n----\nChecking for extra seq IDs in $file1\n----\n";

my $extraprovidercount = 0;
open EXTRA1, ">$outdir/$provider.extra.ids" or die $!;
for my $id (sort keys %providerhash) {
    if (not exists $exportedhash{$id}) {
        print EXTRA1 "$id\n";
        $extraprovidercount++;
    }
}
close EXTRA1;
print LOG "Num of extra      seq IDs in $provider: " . $extraprovidercount . "\n";

####

#print LOG "\n----\nChecking for extra seq IDs in $file2\n----\n";

my $extraexportedcount = 0;
open EXTRA2, ">$outdir/database.extra.ids" or die $!;
for my $id (sort keys %exportedhash) {
    if (not exists $providerhash{$id}) {
        print EXTRA2 "$id\n";
        $extraexportedcount++;
    }
}
close EXTRA2;
print LOG "Num of extra      seq IDs in database: " . $extraexportedcount . "\n";

####

print LOG "Checking sequences with same ID in $provider and database\n";

my @identical;
my @substrproviderofexported;
my @substrexportedofprovider;
my $substrproviderofexportedtostop = 0;
my $substrexportedofprovidertostop = 0;
my @different;

for my $id (sort keys %providerhash) {
    if (exists $exportedhash{$id}) {
        if ($providerhash{$id} eq $exportedhash{$id}) {
            push @identical, $id;
        }
        elsif (index ($providerhash{$id},$exportedhash{$id}) >=0) {
            my $providertostop = $providerhash{$id};
            $providertostop =~ s/\*.*//;
            my $exportedtostop = $exportedhash{$id};
            $exportedtostop =~ s/\*.*//;
            push @substrexportedofprovider, $id . "\t" . distance($providerhash{$id},$exportedhash{$id}) . "\t" . distance($providertostop,$exportedtostop) . "\t" . index ($providerhash{$id},$exportedhash{$id});
            $substrexportedofprovidertostop++ if (index($providertostop,$exportedtostop));
        }
        elsif (index ($exportedhash{$id},$providerhash{$id}) >=0) {
            my $providertostop = $providerhash{$id};
            $providertostop =~ s/\*.*//;
            my $exportedtostop = $exportedhash{$id};
            $exportedtostop =~ s/\*.*//;
            push @substrproviderofexported, $id . "\t" . distance($providerhash{$id},$exportedhash{$id}) . "\t" . distance($providertostop,$exportedtostop) . "\t" . index ($exportedhash{$id},$providerhash{$id});
            $substrproviderofexportedtostop++ if (index($exportedtostop,$providertostop));
        }
        else {
            push @different,  $id . "\t" . distance($providerhash{$id},$exportedhash{$id});
        }
    }
}

print LOG "Identical: " . scalar(@identical) . "\n";
print LOG "Sequence in database is substr of sequence in $provider: $substrexportedofprovidertostop (" . scalar(@substrexportedofprovider) . ")\n";
print LOG "Sequence in $provider is substr of sequence in database: $substrproviderofexportedtostop (" . scalar(@substrproviderofexported) . ")\n";
print LOG "Different: " . scalar(@different) . "\n";

open  OUT, ">$outdir/identical.ids" or die $!;
print OUT  join("\n",@identical)  . "\n" if scalar @identical > 0;
close OUT;
open  OUT, ">$outdir/substrdatabaseofprovider.ids" or die $!;
print OUT  join("\n",@substrexportedofprovider) . "\n" if scalar @substrexportedofprovider > 0;
close OUT;
open  OUT, ">$outdir/substrproviderofdatabase.ids" or die $!;
print OUT  join("\n",@substrproviderofexported) . "\n" if scalar @substrproviderofexported > 0;
close OUT;
open  OUT, ">$outdir/different.ids" or die $!;
print OUT  join("\n",@different) . "\n" if scalar @different > 0;
close OUT;
