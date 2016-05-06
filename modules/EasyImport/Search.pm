#!/usr/bin/perl -w

use strict;
use DBI;

{
  my %strings;
  sub index_seq_regions {
  	my ($core_dbh,$search_dbh,$production_name) = @_;
  	my $sth1 = $core_dbh->prepare("SELECT seq_region_id,name,length FROM seq_region");
  	$sth1->execute();
  	while (my ($sr_id,$name,$length) = $sth1->fetchrow_array()){
  		%strings = ();
  		add_search_term($search_dbh,'seq_region',$sr_id,$name,'seq_region.name',$production_name,$length);
  		my $sth2 = $core_dbh->prepare("SELECT synonym FROM seq_region_synonym WHERE seq_region_id = $sr_id");
  		$sth2->execute();
  		while (my ($synonym) = $sth2->fetchrow_array()){
  			add_search_term($search_dbh,'seq_region',$sr_id,$synonym,'seq_region_synonym.synonym',$production_name);
  		}
  	}
  }

  sub index_genes {
  	my ($core_dbh,$search_dbh,$production_name) = @_;
  	my $sth1 = $core_dbh->prepare("SELECT g.gene_id,g.stable_id,g.description,g.seq_region_start,g.seq_region_end,sr.name,g.display_xref_id FROM gene as g "
  									."JOIN seq_region as sr ON g.seq_region_id = sr.seq_region_id");
  	$sth1->execute();
  	while (my ($g_id,$g_stable_id,$g_desc,$g_start,$g_end,$sr_name,$g_disp_x_id) = $sth1->fetchrow_array()){
  		%strings = ();
  		my $detail = "$sr_name:$g_start-$g_end";
  		add_search_term($search_dbh,'gene',$g_id,$g_stable_id,'gene.stable_id',$production_name,$detail);
  		if ($g_desc){
  			add_search_term($search_dbh,'gene',$g_id,$g_desc,'gene.description',$production_name);
  		}
  		index_xrefs($core_dbh,$search_dbh,$production_name,$g_id,'gene',$g_id,$g_disp_x_id);
  		my $sth2 = $core_dbh->prepare("SELECT transcript_id,stable_id,description,display_xref_id FROM transcript WHERE gene_id = $g_id");
  		$sth2->execute();
  		while (my ($tsc_id,$tsc_stable_id,$tsc_desc,$tsc_disp_x_id) = $sth2->fetchrow_array()){
  			add_search_term($search_dbh,'gene',$g_id,$tsc_stable_id,'transcript.stable_id',$production_name);
  			if ($tsc_desc && (!$g_desc || deslash_and_match($g_desc,$tsc_desc))){
  				add_search_term($search_dbh,'gene',$g_id,$tsc_desc,'transcript.description',$production_name,$tsc_stable_id);
  			}
  			index_xrefs($core_dbh,$search_dbh,$production_name,$g_id,'transcript',$tsc_id,$tsc_disp_x_id);
  			my $sth3 = $core_dbh->prepare("SELECT translation_id, stable_id FROM translation WHERE transcript_id = $tsc_id");
  			$sth3->execute();
  			while (my ($tsl_id,$tsl_stable_id) = $sth3->fetchrow_array()){
  				add_search_term($search_dbh,'gene',$g_id,$tsl_stable_id,'translation.stable_id',$production_name,$tsc_stable_id);
  				index_xrefs($core_dbh,$search_dbh,$production_name,$g_id,'translation',$tsl_id);
  			}

  		}

  	}

  }

  sub index_xrefs {
  	my ($core_dbh,$search_dbh,$production_name,$g_id,$table,$id,$disp_x_id) = @_;
  	my $sth1 = $core_dbh->prepare("SELECT x.xref_id, x.dbprimary_acc,x.display_label,x.description,edb.db_display_name
  									FROM xref x
  									JOIN object_xref o
  									ON x.xref_id = o.xref_id
  									JOIN external_db edb
  									ON x.external_db_id = edb.external_db_id
  									WHERE o.ensembl_id = $id AND o.ensembl_object_type = '$table'");
  	$sth1->execute();
  	while (my ($x_id,$acc,$label,$desc,$edb_name) = $sth1->fetchrow_array()){
  		add_search_term($search_dbh,'gene',$g_id,$acc,$table.'.xref.dbprimary_acc',$production_name,$edb_name);
  		if ($label && $label ne $acc){
  			add_search_term($search_dbh,'gene',$g_id,$label,$table.'.xref.display_label',$production_name,$acc);
  		}
  		if ($disp_x_id && $x_id == $disp_x_id){
  			add_search_term($search_dbh,'gene',$g_id,$label,$table.'.display_name',$production_name);
  		}
  		add_search_term($search_dbh,'gene',$g_id,$desc,$table.'.xref.description',$production_name,$acc) if $desc;
  		my $sth2 = $core_dbh->prepare("SELECT synonym FROM external_synonym WHERE xref_id = $x_id");
  		$sth2->execute();
  		while (my ($synonym) = $sth2->fetchrow_array()){
  			add_search_term($search_dbh,'gene',$g_id,$synonym,$table.'.xref.external_synonym',$production_name,$acc);
  		}
  	}
  }

  sub add_search_term {
  	my ($dbh,$feature_type,$core_id,$string,$string_type,$production_name,$detail) = @_;
  	$string =~ s/\s+$//;
  	if ($feature_type ne 'seq_region'){
  		if (!$strings{$string}){
  			$dbh->do("INSERT INTO multi (feature_type,core_id,string,string_type,production_name)"
  						."VALUES (".$dbh->quote($feature_type)
  						.",".$core_id
  						.",".$dbh->quote($string)
  						.",".$dbh->quote($string_type)
  						.",".$dbh->quote($production_name)
  						.")");
  			my $feature_id = $dbh->last_insert_id(undef,undef,undef,undef);
  			if ($string !~ m/\s/){
  				$dbh->do("INSERT INTO multi_32 (feature_id,feature_type,core_id,string,string_type,production_name)"
  						."VALUES (".$feature_id
  						.",".$dbh->quote($feature_type)
  						.",".$core_id
  						.",".$dbh->quote($string)
  						.",".$dbh->quote($string_type)
  						.",".$dbh->quote($production_name)
  						.")");
  			}
  			else {
  				$dbh->do("INSERT INTO multi_255 (feature_id,feature_type,core_id,string,string_type,production_name)"
  						."VALUES (".$feature_id
  						.",".$dbh->quote($feature_type)
  						.",".$core_id
  						.",".$dbh->quote($string)
  						.",".$dbh->quote($string_type)
  						.",".$dbh->quote($production_name)
  						.")");
  			}
  		}
  	}
  	if ($feature_type ne 'species'){
  	$dbh->do("INSERT INTO $production_name (feature_type,core_id,string,string_type,production_name,detail)"
  						."VALUES (".$dbh->quote($feature_type)
  						.",".$core_id
  						.",".$dbh->quote($string)
  						.",".$dbh->quote($string_type)
  						.",".$dbh->quote($production_name)
  						.",".$dbh->quote($detail)
  						.")");
  		my $feature_id = $dbh->last_insert_id(undef,undef,undef,undef);
  		if (!$strings{$string}){
  		 	if ($string !~ m/\s/){
  				$dbh->do("INSERT INTO $production_name"."_32 (feature_id,feature_type,core_id,string,string_type,production_name)"
  						."VALUES (".$feature_id
  						.",".$dbh->quote($feature_type)
  						.",".$core_id
  						.",".$dbh->quote($string)
  						.",".$dbh->quote($string_type)
  						.",".$dbh->quote($production_name)
  						.")");
  			}
  			else {
  				$dbh->do("INSERT INTO $production_name"."_255 (feature_id,feature_type,core_id,string,string_type,production_name)"
  						."VALUES (".$feature_id
  						.",".$dbh->quote($feature_type)
  						.",".$core_id
  						.",".$dbh->quote($string)
  						.",".$dbh->quote($string_type)
  						.",".$dbh->quote($production_name)
  						.")");
  			}
  		}
  	}
  	$strings{$string} = 1;
  	return;
  }

  sub deslash_and_match {
  	my ($x,$y) = shift;
  	return unless $y;
  	print $y,"\n";
  	$x =~ s/\//\\\//g;
  	$y =~ s/\//\\\//g;
  	return $x =~ m/$y/;
  }

  sub production_name {
  	my $dbh = shift;
  	my $sth = $dbh->prepare("SELECT meta_value FROM meta WHERE meta_key = 'species.production_name'");
  	$sth->execute();
  	if ($sth->rows > 0){
  		return $sth->fetchrow_arrayref()->[0];
  	}
  	return;
  }


  sub setup_search_db {
  	my $dbh = shift;
  	my $params = shift;
    my $dirname = shift;

  	$dbh->do('DROP DATABASE IF EXISTS '.$params->{'DATABASE_SEARCH'}{'NAME'}) || die "ERROR: unable to drop existing [DATABASE_SEARCH] NAME ".$params->{'DATABASE_SEARCH'}{'NAME'}." using provided settings";
  	$dbh->do('CREATE DATABASE '.$params->{'DATABASE_SEARCH'}{'NAME'}) || die "ERROR: unable to create to [DATABASE_SEARCH] NAME ".$params->{'DATABASE_SEARCH'}{'NAME'}." using provided settings";

  	my $file = "$dirname/../sql/search_tables.sql";
  	my $connection_info = 	 '-u'.$params->{'DATABASE_SEARCH'}{'RW_USER'}
  							.' -h'.$params->{'DATABASE_SEARCH'}{'HOST'}
  							.' -P'.$params->{'DATABASE_SEARCH'}{'PORT'}
  							.' -D'.$params->{'DATABASE_SEARCH'}{'NAME'}
  							.' -p'.$params->{'DATABASE_SEARCH'}{'RW_PASS'};
  	system "mysql $connection_info < $file";
  	return $dbh;
  }

  sub search_db_connect {
  	my $params = shift;
    my $dirname = shift;
  	my $dbh = search_db_host_connect($params);
    if ($dbh->do('USE '.$params->{'DATABASE_SEARCH'}{'NAME'})){
      return $dbh;
    }
    $dbh = setup_search_db($dnh,$params,$dirname);
    $dbh->do('USE '.$params->{'DATABASE_SEARCH'}{'NAME'}) || die "ERROR: unable to connect to [DATABASE_SEARCH] NAME ".$params->{'DATABASE_SEARCH'}{'NAME'}." using provided settings";
  	return $dbh;
  }

  sub search_db_host_connect {
  	my $params = shift;
  	my $dsn = "DBI:mysql:host=$params->{'DATABASE_SEARCH'}{'HOST'};port=$params->{'DATABASE_SEARCH'}{'PORT'}";
  	my $dbh = DBI->connect($dsn,"$params->{'DATABASE_SEARCH'}{'RW_USER'}","$params->{'DATABASE_SEARCH'}{'RW_PASS'}") || die "ERROR: unable to connect to [DATABASE_SEARCH] HOST ".$params->{'DATABASE_SEARCH'}{'HOST'}." using provided settings";
  	return $dbh;
  }
}

1;
