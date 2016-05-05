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
use Ensembl_Import;
use Xref_Import;

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
foreach my $subsection (keys %{$params->{'FILES'}}){
	($infiles{$subsection}{'name'},$infiles{$subsection}{'type'}) = fetch_file($params->{'FILES'}{$subsection});
}

## check that an REPEATMASKER file has been specified
die "ERROR: a repeatmasker output file must be specified at [FILES] REPEATMASKER\n" unless $infiles{'REPEATMASKER'};

## connect to the db
my $dbh = core_db_connect($params);

## load the REPEATMASKER file into the database
push @ARGV,$infiles{'REPEATMASKER'}{'name'};
my $hits = read_repeatmasker($dbh,$params,15); # repeatmasker analysis_id should always be 15






__DATA__
   SW   perc perc perc  query           position in query           matching           repeat               position in repeat
score   div. del. ins.  sequence        begin   end        (left)   repeat             class/family       begin  end    (left)     ID

   15   13.7  0.0  0.0  NC_025322.1         282     305   (15709) + (T)n               Simple_repeat           1     24    (0)      1
   11   24.1  0.0  0.0  NC_025322.1        1658    1691   (14323) + (TTATAG)n          Simple_repeat           1     34    (0)      2
   14   15.6  2.6  2.6  NC_025322.1        2578    2615   (13399) + (TATTATC)n         Simple_repeat           1     38    (0)      3
   21   15.9  2.7  8.4  NC_025322.1        3156    3230   (12784) + (TTAATA)n          Simple_repeat           1     71    (0)      4
   16   17.4  4.3  2.1  NC_025322.1        3286    3332   (12682) + (TTACT)n           Simple_repeat           1     48    (0)      5
  252   15.8  0.0  0.0  NW_011953800.1     1757    1794      (54) + rnd-4_family-171   Unknown                31     68  (308) 607217
  936    5.8 14.4  4.0  NW_011953801.1        1     180    (1662) C rnd-5_family-4480  Unknown             (230)    199      2 607218
  292    2.4  0.0  0.0  NW_011953801.1      186     227    (1615) + rnd-4_family-2940  Unknown               222    263    (0) 607219
 3208    5.3  3.5  3.0  NW_011953801.1      231     658    (1184) C rnd-4_family-302   LINE/Jockey         (217)    430      1 607220
  584    5.6  0.0  0.0  NW_011953801.1      660     730    (1112) C rnd-4_family-221   LINE/RTE-RTE        (235)   1006    936 607221
  529   14.1  3.0  2.0  NW_011953801.1      730     830    (1012) + rnd-4_family-1649  Unknown               814    915  (227) 607222 *
 1470   19.3  2.1  1.8  NW_011953801.1      741    1077     (765) C rnd-5_family-761   LINE/RTE-RTE        (146)    360     23 607223
 1220   12.6  6.7  7.5  NW_011953801.1     1342    1624     (218) C rnd-5_family-36    LINE/Jockey         (201)    658    378 607224
 2145    1.5  0.7  3.0  NW_011953802.1        1     270     (909) + rnd-4_family-2543  Unknown               738   1001  (227) 607225
   16    0.0  0.0  0.0  NW_011953802.1      927     943     (236) + (A)n               Simple_repeat           1     17    (0) 607226


__END__
mysql> select distinct repeat_type from repeat_consensus;
+-------------------------+
| repeat_type             |
+-------------------------+
| Dust                    |
| Low complexity regions  |
| LTRs                    |
| RNA repeats             |
| Satellite repeats       |
| Simple repeats          |
| Tandem repeats          |
| Type I Transposons/LINE |
| Type I Transposons/SINE |
| Type II Transposons     |
| Unknown                 |
+-------------------------+
