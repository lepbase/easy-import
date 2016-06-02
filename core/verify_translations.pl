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
  'FILES' =>	{ 	'PROTEIN' => 1
          },
  'META' =>	{ 	'SPECIES.DISPLAY_NAME' => 1
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

my $display_name = $params->{'META'}{'SPECIES.DISPLAY_NAME'};
$display_name =~ s/ /_/g;

my $exported = $display_name.'_-_proteins.fa';
my $provider = $infiles{'PROTEIN'}{'name'};

####

#print STDERR "\n----\nChecking for duplicates in $provider\n----\n";

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
        print STDERR "DUPLICATE $id\n";
        push  @providerdups, $id;
    }
    else {
        $providerhash{$id} = $seqstr;
    }
}

print STDERR "Num of unique     seq IDs in $provider: " . scalar (keys %providerhash) . "\n";
print STDERR "Num of duplicated seq IDs in $provider: " . scalar (@providerdups) . "\n";

if (scalar @providerdups > 0) {
    open  DUP1, ">$provider.dups.ids" or die $!;
    print DUP1  join("\n",@providerdups)  . "\n";
}

####

#print STDERR "\n----\nChecking for duplicates in $exported\n----\n";

my $exportedin = Bio::SeqIO->new(
                            -file   => "<$exported",
                            -format => "FASTA",
                            );

my %exportedhash;
my @exporteddups;

while ( my $seq = $exportedin->next_seq) {
    my $id  = $seq->display_id;
    my $seqstr = $seq->seq;
    $seqstr =~ s/\W*$//;
    if (exists $exportedhash{$id}) {
        print STDERR "DUPLICATE $id\n";
        push  @exporteddups, $id;
    }
    else {
        $exportedhash{$id} = $seqstr;
    }
}

print STDERR "Num of unique     seq IDs in $exported: " . scalar (keys %exportedhash) . "\n";
print STDERR "Num of duplicated seq IDs in $exported: " . scalar (@exporteddups) . "\n";

if (scalar @exporteddups > 0) {
    open  DUP1, ">$exported.dups.ids" or die $!;
    print DUP1  join("\n",@exporteddups)  . "\n";
}

####

#print STDERR "\n----\nChecking for extra seq IDs in $file1\n----\n";

my $extraprovidercount = 0;
open EXTRA1, ">$provider.extra.ids" or die $!;
for my $id (sort keys %providerhash) {
    if (not exists $exportedhash{$id}) {
        print EXTRA1 "$id\n";
        $extraprovidercount++;
    }
}
close EXTRA1;
print STDERR "Num of extra      seq IDs in $provider: " . $extraprovidercount . "\n";

####

#print STDERR "\n----\nChecking for extra seq IDs in $file2\n----\n";

my $extraexportedcount = 0;
open EXTRA2, ">$exported.extra.ids" or die $!;
for my $id (sort keys %exportedhash) {
    if (not exists $providerhash{$id}) {
        print EXTRA2 "$id\n";
        $extraexportedcount++;
    }
}
close EXTRA2;
print STDERR "Num of extra      seq IDs in $exported: " . $extraexportedcount . "\n";

####

print STDERR "Checking sequences with same ID in $provider and $exported\n";

my @identical;
my @substrproviderofexported;
my @substrexportedofprovider;
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
        }
        elsif (index ($exportedhash{$id},$providerhash{$id}) >=0) {
            my $providertostop = $providerhash{$id};
            $providertostop =~ s/\*.*//;
            my $exportedtostop = $exportedhash{$id};
            $exportedtostop =~ s/\*.*//;
            push @substrproviderofexported, $id . "\t" . distance($providerhash{$id},$exportedhash{$id}) . "\t" . distance($providertostop,$exportedtostop) . "\t" . index ($exportedhash{$id},$providerhash{$id});
        }
        else {
            push @different,  $id . "\t" . distance($providerhash{$id},$exportedhash{$id});
        }
    }
}

print STDERR "Identical: " . scalar(@identical) . "\n";
print STDERR "Sequence in $exported is substr of sequence in $provider: " . scalar(@substrexportedofprovider) . "\n";
print STDERR "Sequence in $provider is substr of sequence in $exported: " . scalar(@substrproviderofexported) . "\n";
print STDERR "Different: " . scalar(@different) . "\n";

open  OUT, ">identical.ids" or die $!;
print OUT  join("\n",@identical)  . "\n" if scalar @identical > 0;
open  OUT, ">substrexportedofprovider.ids" or die $!;
print OUT  join("\n",@substrexportedofprovider) . "\n" if scalar @substrexportedofprovider > 0;
open  OUT, ">substrproviderofexported.ids" or die $!;
print OUT  join("\n",@substrproviderofexported) . "\n" if scalar @substrproviderofexported > 0;
open  OUT, ">different.ids" or die $!;
print OUT  join("\n",@different) . "\n" if scalar @different > 0;
