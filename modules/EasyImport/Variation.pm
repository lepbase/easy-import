#!/usr/bin/perl -w

use strict;
use DBI;

# create the variation database from a template
# populate meta table
sub setup_variation_db {
  my $params = shift;
  if (!$params->{'DATABASE_VARIATION'}){
    $params->{'DATABASE_VARIATION'}{'NAME'} = $params->{'DATABASE_CORE'}{'NAME'};
    $params->{'DATABASE_VARIATION'}{'NAME'} =~ s/_core_/_variation_/;
  }
  my $dbh = variation_db_connect($params);
  $dbh->do('DROP DATABASE IF EXISTS '.$params->{'DATABASE_VARIATION'}{'NAME'}) || die "ERROR: unable to drop existing [DATABASE_VARIATION] NAME ".$params->{'DATABASE_VARIATION'}{'NAME'}." using provided settings";
  $dbh->do('CREATE DATABASE '.$params->{'DATABASE_VARIATION'}{'NAME'}) || die "ERROR: unable to create to [DATABASE_VARIATION] NAME ".$params->{'DATABASE_VARIATION'}{'NAME'}." using provided settings";

  my $file = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-variation/sql/table.sql';
  my $connection_info =    '-u'.$params->{'DATABASE_CORE'}{'RW_USER'}
              .' -h'.$params->{'DATABASE_CORE'}{'HOST'}
              .' -P'.$params->{'DATABASE_CORE'}{'PORT'}
              .' -D'.$params->{'DATABASE_VARIATION'}{'NAME'}
              .' -p'.$params->{'DATABASE_CORE'}{'RW_PASS'};
  system "mysql $connection_info < $file";

  $dbh->do('USE '.$params->{'DATABASE_VARIATION'}{'NAME'});
  $dbh->do("INSERT INTO meta (meta_key,meta_value) "
      ."VALUES ('species.production_name',$dbh->quote($params->{'META'}{'SPECIES.PRODUCTION_NAME'}))");
  $dbh->do("INSERT INTO meta (meta_key,meta_value) "
      ."VALUES ('species.division',$dbh->quote($params->{'META'}{'SPECIES.DIVISION'}))");

  return $dbh;
}

sub compara_db_connect {
  my $params = shift;
  my $dsn = "DBI:mysql:host=$params->{'DATABASE_CORE'}{'HOST'};port=$params->{'DATABASE_CORE'}{'PORT'}";
  my $dbh = DBI->connect($dsn,"$params->{'DATABASE_CORE'}{'RW_USER'}","$params->{'DATABASE_CORE'}{'RW_PASS'}") || die "ERROR: unable to connect to [DATABASE_VARIATION] HOST ".$params->{'DATABASE_CORE'}{'HOST'}." using provided settings";
  $dbh->do('CREATE DATABASE IF NOT EXISTS '.$params->{'DATABASE_VARIATION'}{'NAME'}) || die "ERROR: unable to create [DATABASE_VARIATION] NAME ".$params->{'DATABASE_VARIATION'}{'NAME'}." using provided settings";
  $dbh->do('USE '.$params->{'DATABASE_VARIATION'}{'NAME'}) || die "ERROR: unable to connect to [DATABASE_VARIATION] NAME ".$params->{'DATABASE_VARIATION'}{'NAME'}." using provided settings";
  return $dbh;
}

sub core_db_host_connect {
  my $params = shift;
  my $dsn = "DBI:mysql:host=$params->{'DATABASE_CORE'}{'HOST'};port=$params->{'DATABASE_CORE'}{'PORT'}";
  my $dbh = DBI->connect($dsn,"$params->{'DATABASE_CORE'}{'RO_USER'}") || die "ERROR: unable to connect to [DATABASE_CORE] HOST ".$params->{'DATABASE_CORE'}{'HOST'}." using provided settings";
  return $dbh;
}

sub core_db_connect {
  my ($db_name,$params) = @_;
  my $dbh = core_db_host_connect($params);
  $dbh->do('USE '.$db_name) || die "ERROR: unable to connect to ".$db_name." using provided settings";
  return $dbh;
}


1;
