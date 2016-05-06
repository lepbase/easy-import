#!/usr/bin/perl -w

use strict;
use GFFTree;
use DBI;



sub add_analysis {
	# analysis_id - fill in analysis and analysis_description tables
	my ($dbh,$analysis,$provider_name,$provider_url) = @_;
	my $sth = $dbh->prepare("SELECT analysis_id FROM analysis WHERE logic_name LIKE '$analysis'");
    $sth->execute;
    return $sth->fetchrow_arrayref()->[0] if $sth->rows > 0;
  	$dbh->do("INSERT INTO analysis (logic_name) VALUES (".$dbh->quote($analysis).")");
    $sth->execute;
    my $analysis_id = $sth->fetchrow_arrayref()->[0];
	$dbh->do("INSERT INTO analysis_description (analysis_id,description,display_label,web_data) VALUES ($analysis_id,".$dbh->quote("Gene predictions from <a href=\"".$provider_url."\">".$provider_name."</a>").",".$dbh->quote($analysis).",".$dbh->quote('{"caption" => "Genes","colour_key" => "[biotype]","default" => {"MultiBottom" => "collapsed_label","MultiTop" => "gene_label","alignsliceviewbottom" => "as_collapsed_label","contigviewbottom" => "transcript_label","contigviewtop" => "gene_label","cytoview" => "gene_label"},"key" => "'.$analysis_id.'_genes","label_key" => "[biotype]","name" => "Genes"}').")");

    # also add a version of this analysis for ncRNA
    $dbh->do("INSERT INTO analysis (logic_name) VALUES (".$dbh->quote($analysis."_ncRNA").")");
    $sth->execute;
    my $nc_id = $analysis_id + 1;
	$dbh->do("INSERT INTO analysis_description (analysis_id,description,display_label,web_data) VALUES ($nc_id,".$dbh->quote("Gene predictions from <a href=\"".$provider_url."\">".$provider_name."</a>").",".$dbh->quote($analysis).",".$dbh->quote('{"caption" => "ncRNA genes","colour_key" => "[biotype]","default" => {"MultiBottom" => "collapsed_label","MultiTop" => "gene_label","alignsliceviewbottom" => "as_collapsed_label","contigviewbottom" => "transcript_label","contigviewtop" => "gene_label","cytoview" => "gene_label"},"key" => "'.$analysis_id.'_ncgenes","label_key" => "[biotype]","name" => "ncRNA genes"}').")");

    return $analysis_id;
}

sub add_external_db {
	my ($dbh,$dbname) = @_;
	my $sth = $dbh->prepare("SELECT external_db_id FROM external_db WHERE db_name = '$dbname'");
    $sth->execute;
    return $sth->fetchrow_arrayref()->[0] if $sth->rows > 0;
  	$dbh->do("INSERT INTO external_db (db_name) VALUES (".$dbh->quote($dbname).")");
    $sth->execute;
    return $sth->fetchrow_arrayref()->[0];
}


sub get_seq_region_id {
	my ($dbh,$seq_name) = @_;
	my $sth = $dbh->prepare("SELECT seq_region_id FROM seq_region WHERE name = '$seq_name'");
    $sth->execute;
    return $sth->fetchrow_arrayref()->[0] if $sth->rows > 0;

	$sth = $dbh->prepare("SELECT seq_region_id FROM seq_region_synonym WHERE synonym = '$seq_name'");
    $sth->execute;
    return $sth->fetchrow_arrayref()->[0] if $sth->rows > 0;

	return 0;
}

sub add_gene {
	my ($dbh,$gene) = @_;
	my $sth = $dbh->prepare("SELECT gene_id FROM gene WHERE seq_region_id = ".$gene->attributes->{_seq_region_id}." AND seq_region_start = ".$gene->attributes->{_start}." AND seq_region_end = ".$gene->attributes->{_end}." AND seq_region_strand = ".$gene->attributes->{_strand}." AND stable_id LIKE ".$dbh->quote($gene->{attributes}->{stable_id}));
    $sth->execute;
    return $sth->fetchrow_arrayref()->[0] if $sth->rows > 0;
    my $description = $gene->attributes->{description} ? $dbh->quote($gene->attributes->{description}) : 'NULL';
  	$dbh->do("INSERT INTO gene (analysis_id,seq_region_id,seq_region_start,seq_region_end,seq_region_strand,stable_id,source,status,description)
  					 VALUES (".$gene->attributes->{_analysis_id}
  					 		.",".$gene->attributes->{_seq_region_id}
  					 		.",".$gene->attributes->{_start}
  					 		.",".$gene->attributes->{_end}
  					 		.",".$gene->attributes->{_strand}
  					 		.",".$dbh->quote($gene->attributes->{stable_id})
  					 		.",".$dbh->quote($gene->attributes->{_data_source})
  					 		.",".$dbh->quote($gene->attributes->{_status})
  					 		.",".$description
  					 		.")");
    $sth->execute;
    my $gene_id = $sth->fetchrow_arrayref()->[0];
	return $gene_id;
}


sub canonical_transcript {
	my ($dbh,$gene,$sources) = @_;
	my $seq_region_id = $gene->attributes->{_seq_region_id};
	my $analysis_id = $sources->{$gene->attributes->{_source}};
	my ($canonical_transcript_id,$biotype);
	my $trans_stable_id;
	# loop through mrnas to find longest then make/lookup transcript_id
	my @mrna = $gene->by_type('mrna');
	@mrna = $gene->daughters() unless $mrna[0];
	my $longest = 0;
	my $index = -1;
	for (my $i = 0; $i < @mrna; $i++){
		if ($mrna[$i]->_length > $longest){
			$longest = $mrna[$i]->_length;
			$index = $i;
		}
	}
	if ($index >= 0){
		$mrna[$index]->attributes->{_canonical} = 1;
		#$mrna[$index]->attributes->{stable_id} = $mrna[$index]->name() unless $mrna[$index]->attributes->{stable_id};
		($canonical_transcript_id,$biotype) = add_transcript($dbh,$gene,$mrna[$index]);
		return ($canonical_transcript_id,$biotype);
	}
	return (undef,undef);
}

sub add_exon {
	my ($dbh,$gene,$exon) = @_;
	my $sth = $dbh->prepare("SELECT exon_id FROM exon WHERE seq_region_id = ".$gene->attributes->{_seq_region_id}." AND seq_region_start = ".$exon->attributes->{_start}." AND seq_region_end = ".$exon->attributes->{_end}." AND seq_region_strand = ".$exon->attributes->{_strand}." AND (stable_id IS NULL OR stable_id LIKE ".$dbh->quote($exon->mother()->attributes->{stable_id}."%").")");
    $sth->execute;
    return $sth->fetchrow_arrayref()->[0] if $sth->rows > 0;
     $dbh->do("INSERT INTO exon (seq_region_id,seq_region_start,seq_region_end,seq_region_strand,phase,end_phase)
  					 VALUES (".$gene->attributes->{_seq_region_id}
  					 		.",".$exon->attributes->{_start}
  					 		.",".$exon->attributes->{_end}
  					 		.",".$exon->attributes->{_strand}
  					 		.",".$exon->_phase
  					 		.",".$exon->attributes->{_end_phase}
  					 		.")");
    $sth->execute;
    my $exon_id = $sth->fetchrow_arrayref()->[0];
	return $exon_id;

#	mysql> describe exon;
#+-------------------+----------------------+------+-----+---------------------+----------------+
#| Field             | Type                 | Null | Key | Default             | Extra          |
#+-------------------+----------------------+------+-----+---------------------+----------------+
#| exon_id           | int(10) unsigned     | NO   | PRI | NULL                | auto_increment |
#| seq_region_id     | int(10) unsigned     | NO   | MUL | NULL                |                |
#| seq_region_start  | int(10) unsigned     | NO   |     | NULL                |                |
#| seq_region_end    | int(10) unsigned     | NO   |     | NULL                |                |
#| seq_region_strand | tinyint(2)           | NO   |     | NULL                |                |
#| phase             | tinyint(2)           | NO   |     | NULL                |                |
#| end_phase         | tinyint(2)           | NO   |     | NULL                |                |
#| is_current        | tinyint(1)           | NO   |     | 1                   |                |
#| is_constitutive   | tinyint(1)           | NO   |     | 0                   |                |
#| stable_id         | varchar(128)         | YES  | MUL | NULL                |                |
#| version           | smallint(5) unsigned | NO   |     | 1                   |                |
#| created_date      | datetime             | NO   |     | 0000-00-00 00:00:00 |                |
#| modified_date     | datetime             | NO   |     | 0000-00-00 00:00:00 |                |
#+-------------------+----------------------+------+-----+---------------------+----------------+
#13 rows in set (0.00 sec)

}

sub exon_transcript {
	my ($dbh,$exon_id,$transcript_id,$rank) = @_;
	my $sth = $dbh->prepare("SELECT rank FROM exon_transcript WHERE exon_id = ".$exon_id." AND transcript_id = ".$transcript_id);
    $sth->execute;
    return $sth->fetchrow_arrayref()->[0] if $sth->rows > 0;
    $dbh->do("INSERT INTO exon_transcript (exon_id,transcript_id,rank)
  					 VALUES (".$exon_id
  					 		.",".$transcript_id
  					 		.",".$rank
  					 		.")");
    #$sth->execute;

	#return 1;
}

sub add_translation {
	my ($dbh,$mrna,$start_exon,$end_exon,$prot_stable_id) = @_;
	my $transcript_id = $mrna->attributes->{_transcript_id};
	my $start_exon_id = $start_exon->attributes->{exon_id};
	my $end_exon_id = $end_exon->attributes->{exon_id};
	my $seq_start = $start_exon->attributes->{ts_start};
	my $seq_end = $end_exon->attributes->{ts_end};
	my $sth = $dbh->prepare("SELECT translation_id FROM translation WHERE transcript_id = ".$mrna->attributes->{_transcript_id}
																	." AND start_exon_id = ".$start_exon->attributes->{exon_id}
																	." AND end_exon_id = ".$end_exon->attributes->{exon_id}
																	." AND seq_start = ".$seq_start
																	." AND seq_end = ".$seq_end
																	." AND stable_id like ".$dbh->quote($prot_stable_id));
    $sth->execute;
    if ($sth->rows > 0){
    	my $retvals = $sth->fetchrow_hashref();
    	return ($retvals->{'translation_id'});
    }
  	$dbh->do("INSERT INTO translation (transcript_id,seq_start,start_exon_id,seq_end,end_exon_id,stable_id)
  					 VALUES (".$mrna->attributes->{_transcript_id}
  					 		.",".$seq_start
  					 		.",".$start_exon->attributes->{exon_id}
  					 		.",".$seq_end
  					 		.",".$end_exon->attributes->{exon_id}
  					 		.",".$dbh->quote($prot_stable_id)
  					 		.")");
    $sth->execute;
    my $translation_id = $sth->fetchrow_hashref()->{'translation_id'};
	return $translation_id;

#	mysql> describe translation;
#+----------------+----------------------+------+-----+---------------------+----------------+
#| Field          | Type                 | Null | Key | Default             | Extra          |
#+----------------+----------------------+------+-----+---------------------+----------------+
#| translation_id | int(10) unsigned     | NO   | PRI | NULL                | auto_increment |
#| transcript_id  | int(10) unsigned     | NO   | MUL | NULL                |                |
#| seq_start      | int(10)              | NO   |     | NULL                |                |
#| start_exon_id  | int(10) unsigned     | NO   |     | NULL                |                |
#| seq_end        | int(10)              | NO   |     | NULL                |                |
#| end_exon_id    | int(10) unsigned     | NO   |     | NULL                |                |
#| stable_id      | varchar(128)         | YES  | MUL | NULL                |                |
#| version        | smallint(5) unsigned | NO   |     | 1                   |                |
#| created_date   | datetime             | NO   |     | 0000-00-00 00:00:00 |                |
#| modified_date  | datetime             | NO   |     | 0000-00-00 00:00:00 |                |
#+----------------+----------------------+------+-----+---------------------+----------------+
#10 rows in set (0.01 sec)


}

sub add_transcript {
	my ($dbh,$gene,$mrna) = @_;
	my %biotypes = ('mRNA' => 'protein_coding',
					'piRNA' => 'piRNA',
					'tRNA' => 'tRNA',
					'rRNA' => 'rRNA',
					'lncRNA' => 'lncRNA',
					'misc_RNA' => 'misc_RNA');
	my $biotype = $gene->{attributes}->{gene_biotype} || $gene->{attributes}->{biotype} || $biotypes{$mrna->{attributes}->{_type}} || 'misc_RNA';
	$biotype = 'pseudogene' if $gene->{attributes}->{pseudo} && $gene->{attributes}->{pseudo} eq 'true';
	my $sth = $dbh->prepare("SELECT transcript_id FROM transcript WHERE seq_region_id = ".$gene->attributes->{_seq_region_id}." AND seq_region_start = ".$mrna->attributes->{_start}." AND seq_region_end = ".$mrna->attributes->{_end}." AND seq_region_strand = ".$mrna->attributes->{_strand}." AND stable_id LIKE ".$dbh->quote($mrna->{attributes}->{stable_id}));
    $sth->execute;
    if ($sth->rows > 0){
    	return $sth->fetchrow_arrayref()->[0];
  	}
  	my $description = $mrna->attributes->{description} ? $dbh->quote($mrna->attributes->{description}) : 'NULL';
  	$dbh->do("INSERT INTO transcript (gene_id,analysis_id,seq_region_id,seq_region_start,seq_region_end,seq_region_strand,source,biotype,stable_id,status,description)
  					 VALUES (".$gene->attributes->{_gene_id}
  					 		.",".$gene->attributes->{_analysis_id}
  					 		.",".$gene->attributes->{_seq_region_id}
  					 		.",".$mrna->attributes->{_start}
  					 		.",".$mrna->attributes->{_end}
  					 		.",".$mrna->attributes->{_strand}
  					 		.",".$dbh->quote($gene->attributes->{_source})
  					 		.",".$dbh->quote($biotype)
  					 		.",".$dbh->quote($mrna->{attributes}->{stable_id})
  					 		.",".$dbh->quote($gene->attributes->{_status})
  					 		.",".$description
  					 		.")");
    $sth->execute;
    my $transcript_id = $sth->fetchrow_arrayref()->[0];
    return ($transcript_id,$biotype);
}

sub display_xrefs {
	# lookup any gene names/ids in the xrefs table
	# if no matches, insert
	# if matches, insert new names and manage hierarchy
	# also possibly demote current BGI names to xrefs
	# select from xref where display_label like gene name...
	# if there is a match, link to gene id via dependent_xref and object_xref tables
	# otherwise create a new xref entry for this gene
	my ($dbh,$xrefs,$table,$ensembl_id,$analysis_id,$externals,$display_name,$external_synonym) = @_;

	my @xrefs;

	if (ref $xrefs eq 'ARRAY'){
		@xrefs = @$xrefs;
	}
	else {
		$xrefs[0] = $xrefs;
	}

	my $display_xref_id;
	my $xref_id;
	my $name;
	#####################################################
	# TODO: load synonyms differently to multiple xrefs #
	#####################################################

	for (my $i = 0; $i < @xrefs; $i++){
		my $object_xref_id;
		my ($disp,$acc,$external_db_id);
		if ($externals){
			for (my $e = 0; $e < @$externals; $e++){
				my $acc_regex = $externals->[$e]{'acc'};
				if ($xrefs[$i] =~ m/$acc_regex/){
					$acc = $1;
					my $disp_regex = $externals->[$e]{'disp'};
					if ($xrefs[$i] =~ m/$disp_regex/){
						$disp = $1;
					}
					else {
						$disp = $acc;
					}
					$external_db_id = $externals->[$e]{'dbid'};
					last;
				}
			}
		}
		if (!$external_db_id){
			warn "WARNING: No rule to handle xrefs like $xrefs[$i]\n";
			next;
		}

		if ((!$display_name && !$external_synonym) || (!$display_name && $external_synonym && $i == 0) || ($display_name && $xrefs[$i] eq $display_name)){
			# obtain an xref id for this Dbxref
			my $sth_xref = $dbh->prepare("SELECT xref_id FROM xref WHERE dbprimary_acc LIKE ".$dbh->quote($acc)." AND external_db_id = ".$external_db_id);
	   		my $sth_oxref = $dbh->prepare("SELECT object_xref_id FROM object_xref WHERE ensembl_id = ".$ensembl_id." AND ensembl_object_type = ".$dbh->quote($table)." AND analysis_id = ".$analysis_id." AND xref_id = ?");
	    	$sth_xref->execute;
	    	if ($sth_xref->rows > 0){
	    		$xref_id = $sth_xref->fetchrow_arrayref()->[0];

	    		# test to see if this feature already has an object_xref
				$sth_oxref->execute($xref_id);
	    		if ($sth_oxref->rows > 0){
	    			my $arrref = $sth_oxref->fetchrow_arrayref();
	    			$object_xref_id = $arrref->[0];
	    		}
	    	}
	    	else {
	    		$dbh->do('INSERT INTO xref (external_db_id,dbprimary_acc,display_label) VALUES ('
	    			.$external_db_id
	    			.','.$dbh->quote($acc)
	    			.','.$dbh->quote($disp)
	    			.')');
		    	$sth_xref->execute;
		    	$xref_id = $sth_xref->fetchrow_arrayref()->[0];
	    	}
	    	if (!$object_xref_id){
	    		# generate an object_xref_id if one does not already exist
	    		$dbh->do("INSERT INTO object_xref (ensembl_id,ensembl_object_type,xref_id,analysis_id)
	  						 VALUES (".$ensembl_id
	  						 		.",".$dbh->quote($table)
	  						 		.",".$xref_id
	  						 		.",".$analysis_id
	  						 		.")");
	    		$sth_oxref->execute($xref_id);
	    		$object_xref_id = $sth_oxref->fetchrow_arrayref()->[0];
	   		}
	   		$display_xref_id = $xref_id;
	   		$name = $disp;
	   	}
	}
	for (my $i = 0; $i < @xrefs; $i++){
	   	if ($external_synonym){
	   		## insert display_names into external_synonym
	   		my $disp;
	   		if ($externals){
				for (my $e = 0; $e < @$externals; $e++){
					my $acc_regex = $externals->[$e]{'acc'};
					if ($xrefs[$i] =~ m/$acc_regex/){
						my $acc = $1;
						my $disp_regex = $externals->[$e]{'disp'};
						if ($xrefs[$i] =~ m/$disp_regex/){
							$disp = $1;
						}
						else {
							$disp = $acc;
						}
						last;
					}
				}
			}
			if ($disp ne $name){
				my $sth_syn = $dbh->prepare("SELECT xref_id,synonym FROM external_synonym WHERE xref_id = $xref_id AND synonym LIKE ".$dbh->quote($disp));
				$sth_syn->execute;
				if ($sth_syn->rows > 0){
					#
				}
				else {
					$dbh->do("INSERT INTO external_synonym (xref_id,synonym) VALUES ($xref_id,".$dbh->quote($disp).")");
				}
			}
	   	}
    }
	return $display_xref_id;
}

sub identity_xref {
	my ($dbh,$match,$target_name,$transcript) = @_;
	$target_name =~ s/.+://;
	my $external_db_id;
	if ($target_name =~ m/^LOC/){ $external_db_id = 1300; }
	elsif ($target_name =~ m/^NM/){ $external_db_id = 1801; }
	elsif ($target_name =~ m/^XM/){ $external_db_id = 1806; }
	elsif ($target_name =~ m/^NP/){ $external_db_id = 1810; }
	elsif ($target_name =~ m/^XP/){ $external_db_id = 1815; }
	elsif ($target_name =~ m/^NR/){ $external_db_id = 1820; }
	elsif ($target_name =~ m/^XR/){ $external_db_id = 1825; }
	elsif ($target_name =~ m/^NG/){ $external_db_id = 1830; }
	else { die "ERROR: No rule to handle xrefs like ".$target_name."\n" }



	# obtain an xref id for this Dbxref
	my $xref_id;
	my $sth_xref = $dbh->prepare("SELECT max(xref_id) FROM xref WHERE display_label LIKE ".$dbh->quote($target_name)." AND external_db_id = ".$external_db_id);
	$sth_xref->execute;
	if ($sth_xref->rows > 0){
		$xref_id = $sth_xref->fetchrow_arrayref()->[0];
	}
	unless ($xref_id) {
		# need to handle matches to targets that are not in the xref table
		my $dbprimary_acc = $target_name;
		$dbprimary_acc =~ s/\.\d+//;
		$dbh->do("INSERT INTO xref (external_db_id, dbprimary_acc, display_label) VALUES ($external_db_id, ".$dbh->quote($dbprimary_acc).", ".$dbh->quote($target_name).")");
		$sth_xref->execute;
		$xref_id = $sth_xref->fetchrow_arrayref()->[0];
	}
	# need to find the analysis_id
	my $analysis_id;
	my $sth_tsc;
	if ($transcript->{attributes}->{_transcript_id}){
		$sth_tsc = $dbh->prepare("SELECT analysis_id FROM transcript WHERE transcript_id = ".$transcript->{attributes}->{_transcript_id});
	}
	else {
		$sth_tsc = $dbh->prepare("SELECT analysis_id FROM gene WHERE gene_id = ".$transcript->{attributes}->{_gene_id});
	}
	$sth_tsc->execute;
	if ($sth_tsc->rows > 0){
		$analysis_id = $sth_tsc->fetchrow_arrayref()->[0];
	}

	# test to see if this feature already has an object_xref
	my $object_xref_id;
	my $sth_oxref;
	if ($transcript->{attributes}->{_transcript_id}){
		$sth_oxref = $dbh->prepare("SELECT object_xref_id FROM object_xref WHERE ensembl_id = ".$transcript->{attributes}->{_transcript_id}." AND ensembl_object_type LIKE ".$dbh->quote('Transcript')." AND analysis_id = ".$analysis_id." AND xref_id = ".$xref_id);
	}
	else {
		$sth_oxref = $dbh->prepare("SELECT object_xref_id FROM object_xref WHERE ensembl_id = ".$transcript->{attributes}->{_gene_id}." AND ensembl_object_type LIKE ".$dbh->quote('Gene')." AND analysis_id = ".$analysis_id." AND xref_id = ".$xref_id);
	}
	$sth_oxref->execute;
	if ($sth_oxref->rows > 0){
		my $arrref = $sth_oxref->fetchrow_arrayref();
		$object_xref_id = $arrref->[0];
	}
	else {
		print $transcript->{attributes}->{_gene_id},"\n" if $transcript->{attributes}->{_gene_id};
		warn "WARNING: could not find an object_xref_id for $target_name\n";
		return;
	}

	# test to see if this object_xref already has an identity_xref
	my $sth_idxref = $dbh->prepare("SELECT object_xref_id FROM identity_xref WHERE object_xref_id = $object_xref_id");
	$sth_idxref->execute;
	return if $sth_idxref->rows > 0;

	# update the xref info_type to SEQUENCE_MATCH
	simple_update($dbh,'xref',{'info_type' => $dbh->quote('SEQUENCE_MATCH')},{'xref_id' => $xref_id});

	# calculate xref_start, xref_end, ensembl_start, ensembl_end, score and cigar_line
	my ($xref_start, $xref_end, $ensembl_start, $ensembl_end, $score, $cigar_line);
	if ($match->{attributes}->{Target_array}){
		$xref_end = -1;
		$xref_start = 999999999;
		for (my $ i = 0; $i < @{$match->{attributes}->{Target_array}}; $i++){
			$match->{attributes}->{Target_array}->[$i] =~ m/\s(\d+)\s(\d+)/;
			$xref_start = $1 if $1 < $xref_start;
			$xref_end = $2 if $2 > $xref_end;
			if ($match->{attributes}->{_score_array}->[$i]){
				$score += $match->{attributes}->{_score_array}->[$i] if $match->{attributes}->{_score_array}->[$i] =~ m/\d/;
			}
			else {
				$score += $match->{attributes}->{_score} if $match->{attributes}->{_score} =~ m/\d/;
			}
		}
	}
	else {
		$match->{attributes}->{Target} =~ m/\s(\d+)\s(\d+)/;
		$xref_start = $1;
		$xref_end = $2;
		$score = $match->{attributes}->{_score}  if $match->{attributes}->{_score} =~ m/\d/;
	}

	$score = 'NULL' unless $score;


	# fetch the exons for the transcript to help with ensembl_start and ensembl_end
	my (@start_array,@end_array);
	push @start_array,$match->{attributes}->{_start} unless $match->{attributes}->{_start_array};
	push @end_array,$match->{attributes}->{_end} unless $match->{attributes}->{_end_array};
	@start_array = @{$match->{attributes}->{_start_array}} if $match->{attributes}->{_start_array};
	@end_array = @{$match->{attributes}->{_end_array}} if $match->{attributes}->{_end_array};
	for (my $i = 0; $i < @start_array; $i++){
		if ($match->{attributes}->{Gap} && !$match->{attributes}->{Gap_array}){
			$cigar_line .= $match->{attributes}->{Gap};
		}
		elsif ($match->{attributes}->{Gap_array} && $match->{attributes}->{Gap_array}[$i]){
			$cigar_line .= $match->{attributes}->{Gap_array}[$i];
		}

		else {
			$cigar_line .= 'M'.($end_array[$i] - $start_array[$i] + 1);
		}
	}
	$cigar_line =~ s/^M//;
	my @gap_array = split /[\sM]+/,$cigar_line;
	$cigar_line = '';
	my $cigar_m = 0;
	for (my $i = 0; $i <= @gap_array; $i++){
		if ($gap_array[$i] && $gap_array[$i] =~ m/^\d+$/){
			$cigar_m += $gap_array[$i].' ';
		}
		else {
			if ($cigar_m > 0){
				$cigar_line .= 'M'.$cigar_m.' ';
			}
			$cigar_line .= $gap_array[$i].' ' if $gap_array[$i];
			$cigar_m = 0;
		}
	}
	$cigar_line =~ s/\s$//;
	my $length = 0;
	while (my $exon = $transcript->next_feature('exon')){
		if ($exon->{attributes}->{_strand} eq $match->{attributes}->{_strand}){
			if ($exon->{attributes}->{_strand} eq '+1'){
				for (my $i = 0; $i < @start_array; $i++){
					if ($exon->{attributes}->{_start} <= $start_array[$i] && $exon->{attributes}->{_end} >= $end_array[$i]){
						$length += $exon->_length();
						if ($i == 0){
							$ensembl_start = $start_array[$i] - $exon->{attributes}->{_start} + 1;
							$length -= $ensembl_start - 1;
						}
						if ($i + 1 == scalar @start_array){
							if ($exon->{attributes}->{_end} > $end_array[$i]){
								$length -= $exon->{attributes}->{_end} - $end_array[$i];
							}
							$ensembl_end = $ensembl_start + $length - 1;
						}
					}
				}
			}
			else { # TODO edit for reverse strand...
				for (my $i = 0; $i < @start_array; $i++){
					if ($exon->{attributes}->{_end} >= $end_array[$i] && $exon->{attributes}->{_start} <= $start_array[$i]){
						$length += $exon->_length();
						if ($i == 0){
							$ensembl_start = $exon->{attributes}->{_end} - $end_array[$i] + 1;
							$length -= $ensembl_start - 1;
						}
						if ($i + 1 == scalar @start_array){
							if ($exon->{attributes}->{_start} < $start_array[$i]){
								$length -= $start_array[$i] - $exon->{attributes}->{_start};
							}
							$ensembl_end = $ensembl_start + $length - 1;
						}
					}
				}
			}
		}
		else {
			# can't handle this type of alignment yet but it should be very rare
		}
	}

	my ($xref_identity,$ensembl_identity);
	$xref_identity = int($match->{attributes}->{identity} * 100 + 0.5) if $match->{attributes}->{identity};
	$xref_identity = int($match->{attributes}->{pct_identity_ungap} + 0.5) if $match->{attributes}->{pct_identity_ungap} && !$xref_identity;
	$xref_identity = 'NULL' unless $xref_identity;
	$ensembl_identity = int($match->{attributes}->{exon_identity} * 100 + 0.5) if $match->{attributes}->{exon_identity};
	$ensembl_identity = int($match->{attributes}->{pct_identity_ungap} + 0.5) if $match->{attributes}->{pct_identity_ungap} && !$ensembl_identity;
	$ensembl_identity = 'NULL' unless $ensembl_identity;

	# insert values into identity_xref table
	$dbh->do("INSERT INTO identity_xref (object_xref_id,xref_identity,ensembl_identity,xref_start,xref_end,ensembl_start,ensembl_end,cigar_line,score)
  					 VALUES (".$object_xref_id
  					 		.",".$xref_identity
  					 		.",".$ensembl_identity
  					 		.",".$xref_start
  					 		.",".$xref_end
  					 		.",".$ensembl_start
  					 		.",".$ensembl_end
  					 		.",".$dbh->quote($cigar_line)
  					 		.",".$score
  					 		.")");
  	# add an e-value if there is one
  	if ($match->{attributes}->{e_value}){
  		simple_update($dbh,'xref',{'e_value' => $match->{attributes}->{e_value}},{'xref_id' => $xref_id});
  	}

}

sub generate_stable_id {
	my ($dbh,$gene,$prefix) = @_;
	my ($stable_id, $stable_gene_id);
	my $sth = $dbh->prepare("SELECT stable_id FROM gene WHERE gene_id = ".$gene->{attributes}->{_gene_id}." AND stable_id IS NOT NULL");
    $sth->execute;
    if ($sth->rows > 0){
    	$stable_gene_id = $sth->fetchrow_arrayref()->[0];
    	$stable_id = $stable_gene_id;
    	$stable_id =~ s/G0/0/;
    	return $stable_id;
    }
	# make/lookup stable_id
	$sth = $dbh->prepare("SELECT stable_id FROM gene WHERE stable_id LIKE '$prefix"."G%' ORDER BY stable_id DESC LIMIT 1");
    $sth->execute;
    if ($sth->rows > 0){
    	$stable_gene_id = $sth->fetchrow_arrayref()->[0];
		$stable_id = $stable_gene_id;
    	$stable_id =~ s/G0/0/;
    }
    else {
    	$stable_id = $prefix.'00000000';
    }
  	$stable_id =~ m/([1-9][0-9]*)$/;
  	my $prev = $1;
  	$prev = 0 unless $prev;
  	my $next = $prev + 1;
  	$prev = '0'.$prev if length $prev < length $next;
  	$stable_id  =~ s/$prev$/$next/;
  	$stable_gene_id = $stable_id;
	$stable_gene_id =~ s/0/G0/;
	$dbh->do("UPDATE gene SET stable_id = '$stable_gene_id' WHERE gene_id = ".$gene->{attributes}->{_gene_id});

  	return ($stable_id);
}

sub simple_update {
	my ($dbh,$table,$update,$conditions) = @_;
	my $values = "SET ";
	foreach my $field (keys %$update){
		$values .= "$field = $update->{$field}, ";
	}
	$values =~ s/,\s$//;
	my $where = "WHERE ";
	foreach my $field (keys %$conditions){
		$where .= "$field = $conditions->{$field} AND ";
	}
	$where =~ s/\sAND\s$//;
	$dbh->do("UPDATE $table $values $where");

}

sub update_unless {
	my ($dbh,$table,$update,$conditions) = @_;
	my $q = "SELECT * FROM $table WHERE ";
	foreach my $field (keys %$update){
		$q .= "$field IS NOT NULL AND ";
	}
	foreach my $field (keys %$conditions){
		$q .= "$field = $conditions->{$field} AND ";
	}
	$q =~ s/\sAND\s$//;
	my $sth = $dbh->prepare($q);
    $sth->execute;
    if ($sth->rows == 0){
    	simple_update($dbh,$table,$update,$conditions);
    }

}

1;
