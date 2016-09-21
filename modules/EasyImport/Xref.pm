#!/usr/bin/perl -w

use strict;
use DBI;

my %external_db_ids;
my %genes;
my %hits;

sub read_blastp {
	my ($dbh,$params) = @_;
	my $external_db_id = $params->{'XREF'}{'BLASTP'}->[0];
	# TODO make external_db if not exists
	# TODO check that blastp was run with the correct command to give expected column output
	my %transcripts;
	while (<>){
		chomp;
		my ($name,@row) = split /\t/;
		next if $transcripts{$name};
		$transcripts{$name} = 1;
		my $transcript_id = transcript_id($dbh,$name);
		my $acc = my $disp = $row[0];
		$acc =~ s/^sp\|(.*?)\|.*/$1/; # if accession is of form sp|XXX|YYY, extract XXX
		my $desc = $row[13] || undef;
		my $xref_id = xref_id($dbh,$external_db_id,$acc,$disp,$desc,'SEQUENCE_MATCH');
		my $blastp_db = $params->{'XREF'}{'BLASTP'}->[1];
		my $analysis_id = analysis_id($dbh,"BLASTP",$blastp_db); # analysis_id expects [ dbconnection, logic_name, db ]
		my $object_xref_id = object_xref_id($dbh,$xref_id,$transcript_id,$analysis_id);
    object_xref($dbh,$object_xref_id,$xref_id,'IEA') if $external_db_id == 2000;
		my $nident = $row[2] * $row[1] / 100;
		my $xid = $nident / $row[12] * 100;
		my $eid = $nident / $row[11] * 100;
		my $xstart = $row[7] * 3 - 2;
		my $xend = $row[8] * 3 - 2;
		my $estart = $row[5] * 3 - 2;
		my $eend = $row[6] * 3 - 2;
		my $cigar_line = btop2cigar($row[14]);
		my $score = $row[10];
		my $evalue = $row[9];
		identity_xref2($dbh,$object_xref_id,$xid,$eid,$xstart,$xend,$estart,$eend,$cigar_line,$score,$evalue);
		if ($eid > 40 && $params->{'XREF'}{'BLASTP'}[2]){
			my $gene_id = update_description($dbh,$name,$acc,$desc,$params->{'XREF'}{'BLASTP'}[2]);
		}
	}
	return 1;
}

sub read_repeatmasker {
	my ($dbh,$params) = @_;
	my (%repeats,%seq_regions);
  <>;
  <>;
  <>;
	while (<>){
		chomp;
    s/^\s+//;
		my ($score,undef,undef,undef,
        $seq_region_name,$seq_region_start,$seq_region_end,undef,$seq_region_strand,
        $repeat_name,$repeat_class,$repeat_start,$repeat_end,$repeat_left,
        undef,undef) = split /\s+/;
    my $analysis = $repeat_class eq 'Simple_repeat' ? 'Simple_repeat' : 'repeatmasker';
    my $analysis_id = analysis_id($dbh,$analysis,$analysis);
    if ($seq_region_strand eq 'C'){
      $seq_region_strand = -1;
      $repeat_start = $repeat_left;
    }
    else {
      $seq_region_strand = 1;
    }

    my $seq_region_id = $seq_regions{$seq_region_name} || get_seq_region_id($dbh,$seq_region_name);
		$seq_regions{$seq_region_name} = $seq_region_id;
		my $repeat_id = $repeats{$repeat_name} || repeat_consensus($dbh,$repeat_name,$repeat_class);
    $repeats{$repeat_name} = $repeat_id;

    add_repeat_feature($dbh,
                       $seq_region_id,$seq_region_start,$seq_region_end,$seq_region_strand,
                       $repeat_start,$repeat_end,$repeat_id,
                       $analysis_id,
                       $score);

	}
	return 1;
}

sub read_iprscan {
	my ($dbh,$params) = @_;
	my $external_db_id = 1200;
	while (<>){
		chomp;
		my ($name,$hash,$protein_length,$analysis,$hitname,$desc,$start,$end,$evalue,undef,undef,$ipr,undef,undef,$goterms) = split /\t/;
		my ($translation_id,$transcript_id,$gene_id) = translation_id($dbh,$name);
		warn "ERROR: no translation_id for $name\n" and next unless $translation_id;
		my $analysis_id = analysis_id($dbh,$analysis,$analysis); # analysis takes logic_name and db but for interproscan we seem to have
		$ipr ||= ipr($dbh,$hitname);
		if ($ipr){
			$hits{$ipr}{'count'}++;
			$hits{$ipr}{'desc'} = $desc if !$hits{$ipr}{'desc'} || length $hits{$ipr}{'desc'} < length $desc;
			$hits{$ipr}{'genecount'}++ unless $genes{$ipr}{$gene_id};
			$genes{$ipr}{$gene_id}++;
			my $xref_id = xref_id($dbh,$external_db_id,$ipr,$ipr,$desc,'SEQUENCE_MATCH');
			my $object_xref_id = object_xref_id($dbh,$xref_id,$transcript_id,$analysis_id);
		}
		$hitname =~ s/G3DSA://;
		$evalue = $evalue ne '-' ? $evalue : 'NULL';
		protein_feature($dbh,$translation_id,$start,$end,$hitname,$analysis_id,$evalue,$desc);
	}
	return \%hits;
}


sub make_IPtop500_table {
	my ($hits,$params) = @_;
	open HTML,">summary/stats_".ucfirst($params->{'META'}{'SPECIES.PRODUCTION_NAME'})."_IPtop500.html";
print HTML "<table class=\"ss tint fixed data_table no_col_toggle\">
  <colgroup>
    <col width=\"10%\">
    <col width=\"50%\">
    <col width=\"20%\">
    <col width=\"20%\">
  </colgroup>
  <thead>
    <tr><th colspan=\"1\" rowspan=\"1\" class=\"sort_numeric sorting\">No.</th><th colspan=\"1\" rowspan=\"1\" class=\"sort_html sorting\">InterPro name</th><th colspan=\"1\" rowspan=\"1\" class=\"sort_numeric sorting\">Number of Genes</th><th colspan=\"1\" rowspan=\"1\" class=\"sort_position_html sorting\">Number of ensembl hits</th></tr>
  </thead>
  <tbody>
";
	my $bg = 1;
	my $ctr = 1;
	foreach my $ipr (sort { $hits{$b}{'genecount'} <=> $hits{$a}{'genecount'} } keys %hits){

print HTML "    <tr class=\"bg$bg\">
      <td class=\"bold\">$ctr</td>
      <td><a href=\"http://www.ebi.ac.uk/interpro/entry/$ipr\">$ipr</a><br>$hits{$ipr}{'desc'}</td>
      <td><a href=\"/".$params->{'META'}{'SPECIES.PRODUCTION_NAME'}."/Location/Genome?ftype=Domain;id=$ipr\">$hits{$ipr}{'genecount'}</a></td>
      <td>$hits{$ipr}{'count'}</td>
    </tr>
";
		$bg = 1 ? 2 : 1;
		$ctr++;
	}

print HTML "  </tbody>
</table>
";
	return 1;
}

sub identity_xref2 {
	my ($dbh,$object_xref_id,$xid,$eid,$xstart,$xend,$estart,$eend,$cigar_line,$score,$evalue) = @_;

	# test to see if this object_xref already has an identity_xref
	my $sth_idxref = $dbh->prepare("SELECT object_xref_id FROM identity_xref WHERE object_xref_id = $object_xref_id");
	$sth_idxref->execute;
	return if $sth_idxref->rows > 0;

	$dbh->do("INSERT INTO identity_xref (object_xref_id,xref_identity,ensembl_identity,xref_start,xref_end,ensembl_start,ensembl_end,cigar_line,score,evalue)
					 VALUES (".$object_xref_id
					 		.",".$xid
					 		.",".$eid
					 		.",".$xstart
					 		.",".$xend
					 		.",".$estart
					 		.",".$eend
					 		.",".$dbh->quote($cigar_line)
					 		.",".$score
					 		.",".$evalue
					 		.")");

}

sub update_description {
	my ($dbh,$name,$acc,$desc,$source) = @_;

	# test to see if the gene for this transcript already has a description
  my $sth = $dbh->prepare("SELECT g.gene_id,g.description FROM gene AS g JOIN transcript AS t ON g.gene_id = t.gene_id WHERE t.stable_id LIKE ".$dbh->quote($name));
  $sth->execute();
  my $values;

  if ($sth->rows == 0){
      my $sth2 = $dbh->prepare("SELECT tsc.stable_id FROM transcript AS tsc JOIN translation AS tsl ON tsc.transcript_id = tsl.transcript_id WHERE tsl.stable_id LIKE ".$dbh->quote($name));
      $sth2->execute();
      return undef if $sth2->rows == 0;
      my $tsc_stable_id = $sth2->fetchrow_arrayref()->[0];
      my $sth = $dbh->prepare("SELECT g.gene_id,g.description FROM gene AS g JOIN transcript AS t ON g.gene_id = t.gene_id WHERE t.stable_id LIKE ".$dbh->quote($tsc_stable_id));
      $sth->execute();
      return undef if $sth->rows == 0;
      $values = $sth->fetchrow_arrayref();
  } else {
      $values = $sth->fetchrow_arrayref();
  }


  my $gene_id = $values->[0];
  my $description = $values->[1] || 'NULL';

  return $gene_id unless $description eq 'Unknown function' or $description eq 'NULL' or $description =~ m/Source:UniProtKB\/TrEMBL/;

  $desc .= " [Source:$source;Acc:" . $acc . "]";
  $dbh->do("UPDATE gene SET description = ".$dbh->quote($desc)." WHERE gene_id = $gene_id");

  return $gene_id;

}

sub btop2cigar {
	my $btop = shift;
	my @arr = split /(\d+)/,$btop;
	my @exp;
	for (my $i = 0; $i < @arr; $i++){
		if ($arr[$i] =~ m/\d/){
			push @exp,$arr[$i];
		}
		else {
			push @exp,split/(.{2})/,$arr[$i];
		}
	}
	my $cigar = '';
	my ($match,$ins,$del);
	for (my $i = 0; $i < @exp; $i++){
		if ($exp[$i] =~ m/\d/){
			$cigar .= $del.'D' if $del;
			$del = undef;
			$cigar .= $ins.'I' if $ins;
			$ins = undef;
			$match += $exp[$i];
		}
		elsif ($exp[$i] =~ m/\w\w/){
			$cigar .= $del.'D' if $del;
			$del = undef;
			$cigar .= $ins.'I' if $ins;
			$ins = undef;
			$match++;
		}
		elsif ($exp[$i] =~ m/^-/){
			$cigar .= $match.'M' if $match;
			$match = undef;
			$cigar .= $del.'D' if $del;
			$del = undef;
			$ins++;
		}
		elsif ($exp[$i] =~ m/-$/){
			$cigar .= $match.'M' if $match;
			$match = undef;
			$cigar .= $ins.'I' if $ins;
			$ins = undef;
			$del++;
		}
	}
	$cigar .= $match.'M' if $match;
	$cigar .= $del.'D' if $del;
	$cigar .= $ins.'I' if $ins;
	return $cigar;
}

sub transcript_id {
	my ($dbh,$name) = @_;
	my $transcript_id;
	my $sth_transcript = $dbh->prepare("SELECT transcript_id FROM transcript WHERE stable_id LIKE ".$dbh->quote($name));
	$sth_transcript->execute;
	if ($sth_transcript->rows > 0){
		$transcript_id = $sth_transcript->fetchrow_arrayref()->[0];
	}
	else {
		my $sth_translation = $dbh->prepare("SELECT transcript_id FROM translation WHERE stable_id LIKE ".$dbh->quote($name));
		$sth_translation->execute;
		if ($sth_translation->rows > 0){
			$transcript_id = $sth_translation->fetchrow_arrayref()->[0];
		}
	}
	return $transcript_id;
}

sub translation_id {
	my ($dbh,$name) = @_;
	my @array;
	my $sth_translation = $dbh->prepare("SELECT tl.translation_id, tl.transcript_id, tc.gene_id FROM translation AS tl JOIN transcript AS tc ON tl.transcript_id = tc.transcript_id WHERE tl.stable_id LIKE ".$dbh->quote($name));
	$sth_translation->execute;
	if ($sth_translation->rows > 0){
		@array = $sth_translation->fetchrow_array();
	}
	return @array;
}

sub analysis_id {
	my ($dbh,$logic_name,$db) = @_;
	my $analysis_id;
	my $sth = $dbh->prepare("SELECT analysis_id FROM analysis WHERE logic_name LIKE ".$dbh->quote($logic_name)." AND db LIKE ".$dbh->quote($db));
	$sth->execute;
	if ($sth->rows > 0){
		$analysis_id = $sth->fetchrow_arrayref()->[0];
	}
	else {
		$dbh->do("INSERT INTO analysis (logic_name,db)
						VALUES (".$dbh->quote($logic_name)
								.",".$dbh->quote($db)
						 		.")");
		$sth->execute;
		$analysis_id = $sth->fetchrow_arrayref()->[0];
    my %descriptions = (
      'Simple_repeat' => [ 'Simple repeat predictions from RepeatMasker',
              'Simple repeats'
            ],
      'repeatmasker' => [ 'Repeat predictions from RepeatMasker',
              'RepeatMasker'
            ],
      'BLASTP' => [ 'BLASTP hits to the Uniprot database',
              'BLASTP'
            ],
      'Gene3D' => [ 'Gene3D predictions from InterPro Scan',
              'Gene3D',
              "{'type' => 'domain'}"
            ],
      'SMART' => [ 'SMART predictions from InterPro Scan',
              'SMART',
              "{'type' => 'domain'}"
            ],
      'Pfam' => [ 'Pfam predictions from InterPro Scan',
              'Pfam',
              "{'type' => 'domain'}"
            ],
      'ProSiteProfiles' => [ 'ProSiteProfiles predictions from InterPro Scan',
              'ProSiteProfiles',
              "{'type' => 'domain'}"
            ],
      'SUPERFAMILY' => [ 'SUPERFAMILY predictions from InterPro Scan',
              'SUPERFAMILY',
              "{'type' => 'domain'}"
            ],
      'Phobius' => [ 'Phobius predictions from InterPro Scan',
              'Phobius'
            ],
      'TMHMM' => [ 'TMHMM predictions from InterPro Scan',
              'TMHMM'
            ],
      'ProSitePatterns' => [ 'ProSitePatterns predictions from InterPro Scan',
              'ProSitePatterns',
              "{'type' => 'domain'}"
            ],
      'SignalP_EUK' => [ 'SignalP_EUK predictions from InterPro Scan',
              'SignalP_EUK'
            ],
      'PRINTS' => [ 'PRINTS predictions from InterPro Scan',
              'PRINTS',
              "{'type' => 'domain'}"
            ],
      'TIGRFAM' => [ 'TIGRFAM predictions from InterPro Scan',
              'TIGRFAM',
              "{'type' => 'domain'}"
            ]
    );


    my $arr_ref = $descriptions{$logic_name};
		my $web_data = $arr_ref->[2] ? $dbh->quote($arr_ref->[2]) : 'NULL';
		$dbh->do("INSERT INTO analysis_description (analysis_id,description,display_label,displayable,web_data)
						VALUES (".$analysis_id
								.",".$dbh->quote($arr_ref->[0])
								.",".$dbh->quote($arr_ref->[1])
								.",".1
								.",".$web_data
						 		.")");
	}
	return $analysis_id;
}

sub ipr {
	my ($dbh,$name) = @_;
	my $ipr;
	my $sth = $dbh->prepare("SELECT interpro_ac FROM interpro WHERE id LIKE ".$dbh->quote($name));
	$sth->execute;
	if ($sth->rows > 0){
		$ipr = $sth->fetchrow_arrayref()->[0];
	}
	return $ipr;
}

sub xref_id {
	# lookup any names/ids in the xrefs table
	# if no matches, insert a new xref entry

	my ($dbh,$external_db_id,$accession,$display,$description,$info_type) = @_;

	my $xref_id;
	$accession = $dbh->quote($accession);
	$display = $display ? $dbh->quote($display) : $dbh->quote($accession);
	$description = $description ? $dbh->quote($description) : 'NULL';
	my $sth_xref = $dbh->prepare("SELECT xref_id FROM xref WHERE dbprimary_acc LIKE $accession AND display_label LIKE $display AND external_db_id = ".$external_db_id);
	$sth_xref->execute;
	if ($sth_xref->rows > 0){
		$xref_id = $sth_xref->fetchrow_arrayref()->[0];
	}
	else {
		$dbh->do("INSERT INTO xref (external_db_id,dbprimary_acc,display_label,description,info_type)
						VALUES (".$external_db_id
								.",".$accession
						 		.",".$display
						 		.",".$description
						 		.",".$dbh->quote($info_type)
						 		.")");
		$sth_xref->execute;
		$xref_id = $sth_xref->fetchrow_arrayref()->[0];
	}
	return $xref_id;
}

sub object_xref_id {
	my ($dbh,$xref_id,$transcript_id,$analysis_id) = @_;
	my $object_xref_id;
	my $sth_oxref = $dbh->prepare("SELECT object_xref_id FROM object_xref WHERE ensembl_id = $transcript_id AND ensembl_object_type LIKE ".$dbh->quote('Transcript')." AND analysis_id = $analysis_id AND xref_id = $xref_id");
	$sth_oxref->execute;
	if ($sth_oxref->rows > 0){
		$object_xref_id = $sth_oxref->fetchrow_arrayref()->[0];
	}
	else {
		$dbh->do("INSERT INTO object_xref (xref_id,ensembl_id,ensembl_object_type,analysis_id)
						VALUES (".$xref_id
								.",".$transcript_id
	  					 		.",".$dbh->quote('Transcript')
	  					 		.",".$analysis_id
	  					 		.")");
	    $sth_oxref->execute;
	    $object_xref_id = $sth_oxref->fetchrow_arrayref()->[0];
	}
	return $object_xref_id;
}

sub ontology_xref_id {
	my ($dbh,$object_xref_id,$xref_id,$linkage_type) = @_;
	my $sth = $dbh->prepare("SELECT * FROM ontology_xref WHERE object_xref_id = $object_xref_id AND source_xref_id = $xref_id");
	$sth->execute;
	if ($sth->rows > 0){
		return 1;
	}
	$dbh->do("INSERT INTO ontology_xref (object_xref_id,source_xref_id,linkage_type)
					VALUES (".$object_xref_id
							.",".$xref_id
  					 		.",".$dbh->quote($linkage_type)
  					 		.")");
  $sth->execute;
  return 2;
}

sub protein_feature {
	my ($dbh,$translation_id,$start,$end,$hitname,$analysis_id,$evalue,$desc) = @_;

	# test to see if this protein_feature already exists
	my $sth = $dbh->prepare("SELECT protein_feature_id FROM protein_feature WHERE translation_id = $translation_id AND analysis_id = $analysis_id AND seq_start = $start AND seq_end = $end AND hit_name LIKE ".$dbh->quote($hitname));
	$sth->execute;
	return if $sth->rows > 0;

	$dbh->do("INSERT INTO protein_feature (translation_id,seq_start,seq_end,hit_name,analysis_id,evalue,hit_description)
						 VALUES (".$translation_id
					 		.",".$start
					 		.",".$end
					 		.",".$dbh->quote($hitname)
					 		.",".$analysis_id
					 		.",".$evalue
					 		.",".$dbh->quote($desc)
					 		.")");

}

sub repeat_consensus {
	my ($dbh,$repeat_name,$repeat_class) = @_;
  my $id;
	my $sth = $dbh->prepare("SELECT repeat_consensus_id FROM repeat_consensus WHERE repeat_name = '$repeat_name'");
  $sth->execute;
  if ($sth->rows > 0){
	  $id = $sth->fetchrow_arrayref()->[0];
	}
	else {
    my $repeat_consensus = $repeat_name;
    $repeat_consensus =~ s/^\((\w+)\)n/$1/;
    $repeat_consensus = $repeat_consensus =~ m/^[ACGT]+$/ ? $repeat_consensus : 'N';
    my $repeat_type = repeat_type($repeat_class);
		$dbh->do("INSERT INTO repeat_consensus (repeat_name,repeat_class,repeat_type,repeat_consensus)
						VALUES (".$dbh->quote($repeat_name)
								.",".$dbh->quote($repeat_class)
	  					 		.",".$dbh->quote($repeat_type)
	  					 		.",".$dbh->quote($repeat_consensus)
	  					 		.")");
	    $sth->execute;
	    $id = $sth->fetchrow_arrayref()->[0];
	}
	return $id;
}

sub repeat_type {
  my $repeat_class = shift;
  $repeat_class = substr($repeat_class,0,3);
  my %types = ( 'LIN' => 'Type I Transposons/LINE',
                'SIN' => 'Type I Transposons/SINE',
                'DNA' => 'Type II Transposons',
                'Low' => 'Low complexity regions',
                'LTR' => 'LTRs',
                'rRN' => 'RNA repeats',
                'Sim' => 'Simple repeats',
                'Sat' => 'Satellite repeats',
                'Tan' => 'Tandem repeats',
                'Dus' => 'Dust'
              );
   return $types{$repeat_class} if $types{$repeat_class};
   return 'Unknown';
}

sub add_repeat_feature {
	my ($dbh,$seq_region_id,$seq_region_start,$seq_region_end,$seq_region_strand,
      $repeat_start,$repeat_end,$repeat_id,$analysis_id,$score) = @_;
  $dbh->do("INSERT INTO repeat_feature (seq_region_id,seq_region_start,seq_region_end,
                                        seq_region_strand,repeat_start,repeat_end,
                                        repeat_consensus_id,analysis_id,score)
						VALUES (".$seq_region_id
								 .",".$seq_region_start
	  					 	 .",".$seq_region_end
	  					 	 .",".$seq_region_strand
     						 .",".$repeat_start
     	  				 .",".$repeat_end
     	  				 .",".$repeat_id
     	  				 .",".$analysis_id
     	  				 .",".$score
	  					 	 .")");
	return 1;
}


1;
