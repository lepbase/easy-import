#!/usr/bin/perl -w

use strict;
use DBI;
use GFFTree;
use GFF_Import;
use List::Util qw(sum);

sub truncate_seq_tables {
	my ($dbh) = @_;
	my @tables = qw ( coord_system dna genome_statistics seq_region seq_region_attrib seq_region_mapping seq_region_synonym );
	foreach (@tables){
		$dbh->do("TRUNCATE TABLE $_");
	}
	truncate_gene_tables($dbh);
	return 1;
}

sub truncate_gene_tables {
	my ($dbh) = @_;
	my @tables = qw ( associated_xref dependent_xref dna_align_feature exon exon_transcript external_synonym gene gene_archive gene_attrib identity_xref intron_supporting_evidence object_xref ontology_xref peptide_archive protein_align_feature protein_feature stable_id_event supporting_feature transcript transcript_attrib transcript_intron_supporting_evidence transcript_supporting_feature translation translation_attrib xref );
	foreach (@tables){
		$dbh->do("TRUNCATE TABLE $_");
	}
	return 1;
}

sub gff_to_ensembl {
	my ($infile,$dbh,$params) = @_;
	my %sources;
	my $data_source = $params->{'META'}{'GENEBUILD.METHOD'} eq 'import' ? $params->{'META'}{'PROVIDER.NAME'} : 'LepBase';
	my $data_url = $params->{'META'}{'GENEBUILD.METHOD'} eq 'import' ? $params->{'META'}{'PROVIDER.URL'} : 'http://lepbase.org';
	my $total_genes = 0;
	my @externals;
	if ($params->{'DBXREFS'}){
		foreach my $key (keys %{$params->{'DBXREFS'}}){
			my %tmp = ();
			$tmp{'db'} = $key;
			$tmp{'dbname'} = $params->{'DBXREFS'}{$key}[1];
			if ($params->{'DBXREFS'}{$key}[0] == -1){
				# lookup/add a new external_db_id
				$tmp{'dbid'} = add_external_db($dbh,$tmp{'dbname'});
			}
			else {
				$tmp{'dbid'} = $params->{'DBXREFS'}{$key}[0];
			}
			$tmp{'acc'} = $params->{'DBXREFS'}{$key}[2];
			$tmp{'disp'} = $params->{'DBXREFS'}{$key}[3] ? $params->{'DBXREFS'}{$key}[3]: $params->{'DBXREFS'}{$key}[2];
			$tmp{'acc'} =~ s/^\///;
			$tmp{'disp'} =~ s/^\///;
			$tmp{'acc'} =~ s/\/$//;
			$tmp{'disp'} =~ s/\/$//;
			push @externals,\%tmp;
		}
	}
	my @ontologies;
	if ($params->{'ONTOLOGY_TERMS'}){
		foreach my $key (keys %{$params->{'ONTOLOGY_TERMS'}}){
			my %tmp = ();
			$tmp{'db'} = $key;
			$tmp{'dbname'} = $params->{'ONTOLOGY_TERMS'}{$key}[1];
			if ($params->{'ONTOLOGY_TERMS'}{$key}[0] == -1){
				# lookup/add a new external_db_id
				$tmp{'dbid'} = add_external_db($dbh,$tmp{'dbname'});
			}
			else {
				$tmp{'dbid'} = $params->{'ONTOLOGY_TERMS'}{$key}[0];
			}
			$tmp{'acc'} = $params->{'ONTOLOGY_TERMS'}{$key}[2];
			$tmp{'disp'} = $params->{'ONTOLOGY_TERMS'}{$key}[3] ? $params->{'ONTOLOGY_TERMS'}{$key}[3]: $params->{'ONTOLOGY_TERMS'}{$key}[2];
			$tmp{'acc'} =~ s/^\///;
			$tmp{'disp'} =~ s/^\///;
			$tmp{'acc'} =~ s/\/$//;
			$tmp{'disp'} =~ s/\/$//;
			push @ontologies,\%tmp;
		}
	}

	unshift @ARGV,$infile;
	my $gff = GFFTree->new({});
	$gff->name('GFF');
	foreach my $key (keys %{$params->{'GFF'}}){
		my $value = $params->{'GFF'}{$key};
		if (ref $value || ref $value eq 'ARRAY') {
			my @value = @$value;
			my $type = shift @value;
			if ($type eq 'MULTILINE'){
				$gff->multiline(@value);
			}
		}
	}
	while ($gff->parse_chunk('change','region')){
		# modify strands to use +1 or -1
		my @features = $gff->descendants();
		while (my $feature = shift @features){
			$feature->{attributes}->{_strand} = $feature->{attributes}->{_strand}.1 * 1 if $feature->attributes->{_strand}
		}
		# add analyses for any new sources
		my @sources = $gff->by_attribute('_source');
		while (my $source = shift @sources){
			if (!$sources{$source->attribute->{_source}}){
				$sources{$source->attribute->{_source}} = add_analysis($dbh,$source->attribute->{_source},$data_source,$data_url);
			}
		}
		my @gene = $gff->by_type('gene');
		for (my $g = 0; $g < @gene; $g++){

			my $analysis_id = $sources{$gene[$g]->attributes->{_source}};
			$gene[$g]->attributes->{_analysis_id} = $analysis_id;

			my $seq_region_id = get_seq_region_id($dbh,$gene[$g]->attributes->{_seq_name});
			if (!$seq_region_id){
				## skip chunk if sequence cannot be found in the database
				warn "WARNING: no seq_region in the database matches '".$gene[$g]->attributes->{_seq_name}."'\n";
				last;
			}
			$gene[$g]->attributes->{_seq_region_id} = $seq_region_id;

			my $seq_region_start = $gene[$g]->attributes->{_start};
			my $seq_region_end = $gene[$g]->attributes->{_end};
			my $seq_region_strand = $gene[$g]->attributes->{_strand};

			$gene[$g]->attributes->{_data_url} = $data_url;
			$gene[$g]->attributes->{_data_source} = $data_source;
			# TODO: allow for status to be something other than 'PREDICTED';
			$gene[$g]->attributes->{_status} = 'PREDICTED';

			my $gene_id = add_gene($dbh,$gene[$g]);
			$gene[$g]->attributes->{_gene_id} = $gene_id;

			# load dbxrefs
			if ($gene[$g]->{attributes}->{Dbxref}){
				display_xrefs($dbh,$gene[$g]->attributes->{Dbxref},'gene',$gene_id,$analysis_id,\@externals);
			}

			# load ontology_terms
			if ($gene[$g]->{attributes}->{Ontology_term}){
				display_xrefs($dbh,$gene[$g]->attributes->{Ontology_term},'gene',$gene_id,$analysis_id,\@ontologies);
			}


			# load synonyms and display_name
			if ($gene[$g]->{attributes}->{synonym}){
				if ($gene[$g]->{attributes}->{display_name}){
					if (my $display_xref_id = display_xrefs($dbh,$gene[$g]->attributes->{synonym},'gene',$gene_id,$analysis_id,\@externals,$gene[$g]->{attributes}->{display_name},1)){
						simple_update($dbh,'gene',{'display_xref_id' => $display_xref_id},{'gene_id' => $gene_id});
					}
				}
				else {
					display_xrefs($dbh,$gene[$g]->attributes->{synonym},'gene',$gene_id,$analysis_id,\@externals,undef,1);
				}
			}

			# canonical_transcript_id - find out longest mrna then add to transcript table to generate an id
			my ($canonical_transcript_id,$biotype) = canonical_transcript($dbh,$gene[$g],\%sources);

			if ($canonical_transcript_id){
				simple_update($dbh,'gene',{'canonical_transcript_id' => $canonical_transcript_id, 'biotype' => $dbh->quote($biotype)},{'gene_id' => $gene_id});

				# loop through mrnas, adding transcript if not canonical and storing exons and translations
				my $type = 'coding';
				my @mrna = $gene[$g]->by_type('mrna');
				if (!$mrna[0]){
					$type = 'non-coding';
					$gene[$g]->{attributes}->{_analysis_id} = $analysis_id + 1;
					simple_update($dbh,'gene',{'analysis_id' => $gene[$g]->{attributes}->{_analysis_id}},{'gene_id' => $gene_id});
					@mrna = $gene[$g]->daughters();
				}
				my $m = 0;
				while (my $mrna = shift @mrna){
					$m++;
					($mrna->attributes->{_transcript_id},$biotype) = add_transcript($dbh,$gene[$g],$mrna);

					# load dbxrefs
					if ($mrna->{attributes}->{Dbxref}){
						display_xrefs($dbh,$mrna->attributes->{Dbxref},'transcript',$mrna->attributes->{_transcript_id},$analysis_id,\@externals)
					}

					# load ontology_terms
					if ($mrna->{attributes}->{Ontology_term}){
						display_xrefs($dbh,$mrna->attributes->{Ontology_term},'transcript',$mrna->attributes->{_transcript_id},$analysis_id,\@ontologies)
					}

					# load synonyms and display_name
					if ($mrna->{attributes}->{synonym}){
						if ($mrna->{attributes}->{display_name}){
							if (my $mrna_display_xref_id = display_xrefs($dbh,$mrna->attributes->{synonym},'transcript',$mrna->attributes->{_transcript_id},$analysis_id,\@externals,$mrna->{attributes}->{display_name},1)){
								simple_update($dbh,'transcript',{'display_xref_id' => $mrna_display_xref_id},{'transcript_id' => $mrna->attributes->{_transcript_id}});
							}
						}
						else {
							display_xrefs($dbh,$mrna->attributes->{synonym},'transcript',$mrna->attributes->{_transcript_id},$analysis_id,\@externals,undef,1);
						}
					}

					my ($cds,@startarr,@ends,@phases);
					if ($type eq 'coding'){
						# loop through cds, to work out phase and end_phase for exons plus start and stop for translation
						$cds = $mrna->by_type('cds');
						if ($cds && $cds->attributes->{_start_array}){
							@startarr = @{$cds->attributes->{_start_array}};
							@ends = @{$cds->attributes->{_end_array}};
							@phases = @{$cds->attributes->{_phase_array}};
						}
						elsif ($cds) {
							$startarr[0] = $cds->attributes->{_start};
							$ends[0] = $cds->attributes->{_end};
							$phases[0] = $cds->attributes->{_phase};

						}
						else {
							@startarr = ();
						}
					}
					# loop through exons,
					# to work out phase and end_phase for exons plus start and stop for translation
					# add to database and store exon_ids
					my @exon = $mrna->by_type('exon');
					my ($first,$last) = (-1,-1);
					my @starts;
					my @adjust = (0,1,2);
					if ($params->{'MODIFY'}{'INVERT_PHASE'}){
						@adjust = (0,2,1);
					}
					for (my $p = 0; $p < @phases; $p++){
						$phases[$p] = $adjust[$phases[$p]];
					}
					my $prot_stable_id = $mrna->attributes->{translation_stable_id} if $cds;
					for (my $e = 0; $e < @exon; $e++){
						$exon[$e]->attributes->{_phase} = -1;
						$exon[$e]->attributes->{_end_phase} = -1;
						for (my $c = 0; $c < @startarr; $c++){
							if ($exon[$e]->attributes->{_start} <= $startarr[$c] && $exon[$e]->attributes->{_end} >= $ends[$c]){
								# exon contains the cds so use cds phase information to set exon phase
								if ($exon[$e]->attributes->{_start} == $startarr[$c]){
									if ($exon[$e]->attributes->{_strand} > 0){
										$exon[$e]->_phase($phases[$c]);
										$exon[$e]->attributes->{ts_start} = 1;
									}
									else {
										$exon[$e]->attributes->{_end_phase} = ($phases[$c] + $ends[$c] - $startarr[$c] + 1) % 3;
										$exon[$e]->attributes->{ts_end} = $exon[$e]->_length();
									}
								}
								elsif ($exon[$e]->attributes->{_strand} > 0){
									$exon[$e]->attributes->{ts_start} = $startarr[$c] - $exon[$e]->attributes->{_start} + 1;
								}
								elsif ($exon[$e]->attributes->{_strand} < 0){
									$exon[$e]->attributes->{ts_end} = $exon[$e]->attributes->{_end} - $startarr[$c] + 1;
								}
								if ($exon[$e]->attributes->{_end} == $ends[$c]){
									if ($exon[$e]->attributes->{_strand} > 0){
										$exon[$e]->attributes->{_end_phase} = ($phases[$c] + $ends[$c] - $startarr[$c] + 1) % 3;
										$exon[$e]->attributes->{ts_end} = $exon[$e]->_length();
									}
									else {
										$exon[$e]->_phase($phases[$c]);
										$exon[$e]->attributes->{ts_start} = 1;
									}
								}
								elsif ($exon[$e]->attributes->{_strand} > 0){
									$exon[$e]->attributes->{ts_end} = $ends[$c] - $exon[$e]->attributes->{_start} + 1;
								}
								elsif ($exon[$e]->attributes->{_strand} < 0){
									$exon[$e]->attributes->{ts_start} = $exon[$e]->attributes->{_end} - $ends[$c] + 1;
								}
								$first = $e+1 unless $first < $e+1;
								$first = $e+1 unless $first > -1;
								$last = $e+1 unless $last > $e+1;

							}
						}
						$exon[$e]->{attributes}->{_transcript_id} = $mrna->attributes->{_transcript_id};
						my $exon_id = add_exon($dbh,$gene[$g],$exon[$e]);
						$exon[$e]->{attributes}->{exon_id} = $exon_id;
						$starts[$e] = $exon[$e]->attributes->{_start};

					}
					my @order = sort { $starts[$a] <=> $starts[$b] } 0 .. $#starts;
					my ($start_exon,$end_exon);
					for (my $r = 0; $r < @order; $r++){
						my $rank = $r+1;
						$rank = @order - $r if $gene[$g]->attributes->{_strand} == -1;
						exon_transcript($dbh,$exon[$order[$r]]->{attributes}->{exon_id},$mrna->attributes->{_transcript_id},$rank);
						update_unless($dbh,'exon',{'stable_id' => $dbh->quote($mrna->attributes->{stable_id}."-E".$rank)},{'exon_id' => $exon[$order[$r]]->{attributes}->{exon_id}});
						if ($first > -1 && $order[$r] == ($first-1)){
							$end_exon = $first;
							$start_exon = $first unless $start_exon;
						}
						if ($last > -1 && $order[$r] == ($last-1)){
							$end_exon = $last;
							$start_exon = "$last" unless $start_exon;
						}
					}
					# then do translation if appropriate
					if ($cds && $first > -1){
						my ($translation_id,$canonical);
						if ($gene[$g]->attributes->{_strand} == 1){
							$translation_id = add_translation($dbh,$mrna,$exon[($start_exon-1)],$exon[($end_exon-1)],$prot_stable_id);
						}
						else {
							$translation_id = add_translation($dbh,$mrna,$exon[($end_exon-1)],$exon[($start_exon-1)],$prot_stable_id);
						}
						simple_update($dbh,'transcript',{'canonical_translation_id' => $translation_id},{'transcript_id' => $mrna->attributes->{_transcript_id}});

					}

				}
				# TODO: incorporate test for pseudogene at some point and add this to the gff
				if ($gene[$g]->attributes->{pseudo} && $gene[$g]->attributes->{pseudo} =~ m/true/i){
					simple_update($dbh,'gene',{'biotype' => $dbh->quote('pseudogene')},{'gene_id' => $gene_id});
				}

			}
		}
	}
	return $total_genes;
}

sub count_rows {
	my ($dbh,@rows) = @_;
	foreach (@rows){
		my $sth = $dbh->prepare("SELECT count(*) FROM ".$_);
		$sth->execute;
		if ($sth->rows > 0){
			my $count = $sth->fetchrow_arrayref()->[0];
			print "$count ".$_."s in database\n";
		}
	}
	return 1;
}

sub rewrite_gff {
	my ($params,$infile,$properties,$stable_ids_ref,$desc_ref,$names_ref,$tr_properties,$tr_stable_ids_ref,$tr_desc_ref,$tr_names_ref,$tl_stable_ids_ref) = @_;
	my (%ids,%dups);
	my ($stable_id_location,$stable_id_regex,$stable_id_substitution) = @$stable_ids_ref;
	my ($index,$desc_location,$desc_regex,$desc_substitution) = @$desc_ref if $desc_ref;
	my ($name_index,$name_location,$name_regex,$name_substitution) = @$names_ref if $names_ref;
	my ($tr_stable_id_location,$tr_stable_id_regex,$tr_stable_id_substitution) = @$tr_stable_ids_ref if $tr_stable_ids_ref;
	my ($tr_index,$tr_desc_location,$tr_desc_regex,$tr_desc_substitution) = @$tr_desc_ref if $tr_desc_ref;
	my ($tr_name_index,$tr_name_location,$tr_name_regex,$tr_name_substitution) = @$tr_names_ref if $tr_names_ref;
	my ($tl_stable_id_location,$tl_stable_id_regex,$tl_stable_id_substitution) = @$tl_stable_ids_ref if $tl_stable_ids_ref;
	# create new gff tree object
	# add lots of validation conditions
	my $filename = $infile->{'name'};
	if ($params->{'GFF'}{'SORT'} && $filename !~ m/.sorted/){
		$filename .= ".sorted";
		if (!-e $filename){
			system "cat ".$infile->{'name'}." | sort > ".$filename.".tmp";
			system "mv ".$filename.".tmp ".$filename;
		}
	}
	unshift @ARGV,$filename;
	open OUT, ">".$filename.".gff.tmp";
	open EXC, ">".$filename.".exception.gff.tmp";
	my $gff = GFFTree->new({});
	$gff->name('GFF');
	$gff->expect_columns(9,'skip');
	$gff->lacks_id('all','make');
	$gff->lacks_id('cds','Name');
  	#$gff->multiline('CDS');
	#$gff->add_expectation('cds','hasSister','exon','make');
	#$gff->add_expectation('mrna','hasParent','gene','force');
	#$gff->add_expectation('trna','hasParent','gene','force');
	#$gff->add_expectation('transcript','hasParent','gene','force');
	#$gff->add_expectation('exon','hasParent','mrna|transcript','force');
	#$gff->add_expectation('cds','hasParent','mrna|transcript','force');
	#$gff->add_expectation('cds|exon|mrna|trna|transcript|gene','<=[_start,_end]','SELF','warn');
	foreach my $key (keys %{$params->{'GFF'}}){
		my $value = $params->{'GFF'}{$key};
		if (ref $value || ref $value eq 'ARRAY') {
			my @value = @$value;
			my $type = shift @value;
			if ($type eq 'MULTILINE'){
				$gff->multiline(@value);
			}
			elsif ($type eq 'LACKS_ID'){
				$gff->lacks_id(@value);
			}
			elsif ($type eq 'EXPECTATION'){
				$gff->add_expectation(@value);
			}
			elsif ($type eq 'MAP_TYPES'){
				$gff->map_types({$value[0] => $value[1]});
			}
			elsif ($key !~ m/^(:?SORT|SPLIT|CHUNK)$/){
				warn "WARNING: No pattern to handle $value in [GFF] $key\n";
			}
		}
	}


	#while ($gff->parse_chunk('separator','###')){
	while ($gff->parse_chunk('change','region')){
		# PRIORITY TODO: skip features by source (if set in the ini)
		$gff->validate_all();
		my @genes = $gff->by_type('gene');
		my @exceptions;
		my @valid;
		while (my $gene = shift @genes){
			my $exception = 0;
			my ($stable_id,$description,$name);
			if ($stable_id_location =~ m/gene->(.+)/){
				$stable_id = $gene->{attributes}->{$1};
			}
			elsif ($stable_id_location =~ m/DAUGHTER->(.+)/){
				my @tmp = $gene->daughters();
				$stable_id = $tmp[0]->{attributes}->{$1};
			}
			elsif ($stable_id_location =~ m/(.+?)->(.+)/){
				my $tmp = $gene->by_type($1);
				$stable_id = $tmp->{attributes}->{$2};
			}
			else {
				warn "WARNING: could not recognise [GENE_STABLE_IDS] $stable_id_location as a location in a file of type ".$infile->{'type'}."\n";
				push @exceptions,$gene;
				next;
			}
			if (!$stable_id){
				warn "WARNING: could not find a stable_id at [GENE_STABLE_IDS] $stable_id_location for gene ".$gene->id."\n";
				push @exceptions,$gene;
				next;
			}

			if ($desc_ref){
				if ($desc_location =~ m/gene->(.+)/){
					$description = $gene->{attributes}->{$1};
				}
 			    elsif ($desc_location =~ m/DAUGHTER->(.+)/){
				    my @tmp = $gene->daughters();
				    $description = $tmp[0]->{attributes}->{$1};
			    }
				elsif ($desc_location =~ m/(.+?)->(.+)/){
					my $tmp = $gene->by_type($1);
					$description = $tmp->{attributes}->{$2};
				}
				else {
					die "ERROR: could not recognise [GENE_DESCRIPTIONS] $desc_location as a location in a file of type ".$infile->{'type'}."\n";
				}
				($stable_id,$description) = _match_property_to_stable_id('description',$properties,$index,$stable_id,$description,$stable_id_regex,$stable_id_substitution,$desc_regex,$desc_substitution);
				if ($ids{'gene'}{$stable_id}){
					$dups{'gene'}{$stable_id}++;
				}
				$ids{'gene'}{$stable_id}++;
				if ($description){
					$gene->{attributes}->{'description'} = $description;
				}
			}
			elsif ($properties->{$stable_id}{'description'}{'values'}[0]){
				$gene->{attributes}->{'description'} = $properties->{$stable_id}{'description'}{'values'}[0];
			}

			if ($names_ref){
				if ($name_location =~ m/gene->(.+)/){
					$name = $gene->{attributes}->{$1};
				}
				elsif ($name_location =~ m/(.+?)->(.+)/){
					my $tmp = $gene->by_type($1);
					$name = $tmp->{attributes}->{$2};
				}
				else {
					die "ERROR: could not recognise [GENE_DESCRIPTIONS] $name_location as a location in a file of type ".$infile->{'type'}."\n";
				}
				## test if $name is actually a list here and loop through all;
				if ($name){
					if (! ref $name || ref $name ne 'ARRAY') {
  						my @tmp = split /,/,$name;
						$name = \@tmp;
					}
					my ($new_stable_id,$new_name);
					for (my $i = 0; $i < @$name; $i++){
						($new_stable_id,$new_name) = _match_property_to_stable_id('synonym',$properties,$index,$stable_id,$name->[$i],$stable_id_regex,$stable_id_substitution,$name_regex,$name_substitution);
					}
					$stable_id = $new_stable_id;
				}
			}
			if ($properties->{$stable_id}{'synonym'}{'values'}[0]){
				my $names = $properties->{$stable_id}{'synonym'}{'values'};
				if (!ref $names || ref $names ne 'ARRAY') {
					my @tmp = ($names);
					$names = \@tmp;
				}
				my $synonyms;
				if ($gene->{attributes}->{'synonym'}){
					$synonyms = $gene->{attributes}->{'synonym'};
					if (! ref $synonyms || ref $synonyms ne 'ARRAY') {
  						my @tmp = ($synonyms);
						$synonyms = \@tmp;
					}
				}
				if ($names && $synonyms){
					$gene->{attributes}->{'synonym'} = [@$names,@$synonyms];
				}
				elsif ($names){
					$gene->{attributes}->{'synonym'} = $names;
				}
			}
			if ($properties->{$stable_id}{'synonym'}{'index'} && $properties->{$stable_id}{'synonym'}{'index'} <= 1){
				$gene->{attributes}->{'display_name'} = $gene->{attributes}->{'synonym'}->[0];
			}
			$gene->{attributes}->{'stable_id'} = $stable_id;

			# also do transcript stable_ids and descriptions
			my @tmp = $gene->daughters();
			while (my $transcript = shift @tmp){
				my ($tr_stable_id,$tr_description,$tr_name,$tl_stable_id);
				if ($tr_stable_id_location){
					$tr_stable_id_location =~ m/(.+?)->(.+)/;
					my $type = $1;
					my $attr = $2;
					if ($type eq 'SELF' || $transcript->{attributes}->{_type} =~ m/($type)/){
						$tr_stable_id = $transcript->{attributes}->{$attr};
					}
					elsif ($type =~ m/gene/){
						$tr_stable_id = $gene->{attributes}->{$attr};
					}
					elsif ($tr_stable_id_location =~ m/(.+?)->(.+)/){
						my $tmp = $transcript->by_type($1);
						$tr_stable_id = $tmp->{attributes}->{$2};
					}
					if ($tr_stable_id){
						$tr_stable_id = _match_and_replace($tr_stable_id,$tr_stable_id_regex,$tr_stable_id_substitution);
						if (!$tr_stable_id){
							# no transcript_stable_id could be found so place gene in exceptions file
							warn "WARNING: could not find a transcript stable_id at [TRANSCRIPT_STABLE_IDS] $tr_stable_id_location for gene ".$gene->id." transcript ".$transcript->id."\n";
							#push @exceptions,$gene;
							$exception = 1;
							last;
						}
						if ($ids{'transcript'}{$tr_stable_id}){
							$dups{'transcript'}{$tr_stable_id}++;
						}
						$ids{'transcript'}{$tr_stable_id}++;

						# look for a translation_stable_id
						if ($tl_stable_id_location){
							$tl_stable_id_location =~ m/(.+?)->(.+)/;
							my $type = $1;
							my $attr = $2;
							if ($type eq 'SELF' || $transcript->{attributes}->{_type} =~ m/($type)/){
								$tl_stable_id = $transcript->{attributes}->{$attr};
							}
							elsif ($type =~ m/gene/){
								$tl_stable_id = $gene->{attributes}->{$attr};
							}
							elsif ($tl_stable_id_location =~ m/(.+?)->(.+)/){
								my $tmp = $transcript->by_type($1);
								$tl_stable_id = $tmp->{attributes}->{$2};
							}
							if ($tl_stable_id){
								$tl_stable_id = _match_and_replace($tl_stable_id,$tl_stable_id_regex,$tl_stable_id_substitution);
								if ($ids{'translation'}{$tl_stable_id}){
									$dups{'translation'}{$tl_stable_id}++;
								}
								$ids{'translation'}{$tl_stable_id}++;
							}
						}
					}
					else {
						# no transcript_stable_id could be found so place gene in exceptions file
						warn "WARNING: could not find a transcript stable_id at [TRANSCRIPT_STABLE_IDS] $tr_stable_id_location for gene ".$gene->id." transcript ".$transcript->id."\n";
						#push @exceptions,$gene;
						$exception = 1;
						last;
					}
				}
				else {
					die "ERROR: could not recognise [TRANSCRIPT_STABLE_IDS] $tr_stable_id_location as a location in a file of type ".$infile->{'type'}."\n";
				}
				if ($tr_desc_ref){
					if ($tr_desc_location){
						$tr_desc_location =~ m/(.+?)->(.+)/;
						my $type = $1;
						my $attr = $2;
						if ($type eq 'SELF' || $transcript->{attributes}->{_type} =~ m/($type)/){
							$tr_description = $transcript->{attributes}->{$attr};
						}
						elsif ($type =~ m/gene/){
							$tr_description = $gene->{attributes}->{$attr};
						}
						elsif ($tr_desc_location =~ m/(.+?)->(.+)/){
							my $tmp = $transcript->by_type($1);
							$tr_description = $tmp->{attributes}->{$2};
						}
					}

					if ($tr_desc_location =~ m/gene->(.+)/){
						$tr_description = $gene->{attributes}->{$1};
					}
					elsif ($tr_desc_location =~ m/(.+?)->(.+)/){
						my $tmp = $gene->by_type($1);
						$tr_description = $tmp->{attributes}->{$2};
					}
					($tr_stable_id,$tr_description) = _match_property_to_stable_id('description',$tr_properties,$tr_index,$tr_stable_id,$tr_description,$tr_stable_id_regex,$tr_stable_id_substitution,$tr_desc_regex,$tr_desc_substitution);
					if ($tr_description){
						$transcript->{attributes}->{'description'} = $tr_description;
					}
				}
				elsif ($tr_properties->{$tr_stable_id}{'description'}{'values'}[0]){
					$transcript->{attributes}->{'description'} = $tr_properties->{$tr_stable_id}{'description'}{'values'}[0];
				}
				if ($tr_names_ref){
					if ($tr_name_location){
						$tr_name_location =~ m/(.+?)->(.+)/;
						my $type = $1;
						my $attr = $2;
						if ($type eq 'SELF' || $transcript->{attributes}->{_type} =~ m/($type)/){
							$tr_name = $transcript->{attributes}->{$attr};
						}
						elsif ($type =~ m/gene/){
							$tr_name = $gene->{attributes}->{$attr};
						}
						elsif ($tr_name_location =~ m/(.+?)->(.+)/){
							my $tmp = $transcript->by_type($1);
							$tr_name = $tmp->{attributes}->{$2};
						}
					}
					## test if $tr_name is actually a list here and loop through all;
					if ($tr_name){
						if (! ref $tr_name || ref $tr_name ne 'ARRAY') {
  							my @tmp = split /,/,$tr_name;
							$tr_name = \@tmp;
						}
						my ($new_tr_stable_id,$new_tr_name);
						for (my $i = 0; $i < @$tr_name; $i++){
							($new_tr_stable_id,$new_tr_name) = _match_property_to_stable_id('synonym',$tr_properties,$tr_index,$tr_stable_id,$tr_name->[0],$tr_stable_id_regex,$tr_stable_id_substitution,$tr_name_regex,$tr_name_substitution);
						}
						$tr_stable_id = $new_tr_stable_id;
					}
				}
				if ($tr_properties->{$tr_stable_id}{'synonym'}{'values'}[0]){
					my $names = $tr_properties->{$tr_stable_id}{'synonym'}{'values'};
					if (!ref $names || ref $names ne 'ARRAY') {
						my @tmp = ($names);
						$names = \@tmp;
					}
					my $synonyms;
					if ($transcript->{attributes}->{'synonym'}){
						$synonyms = $transcript->{attributes}->{'synonym'};
						if (! ref $synonyms || ref $synonyms ne 'ARRAY') {
  							my @tmp = ($synonyms);
							$synonyms = \@tmp;
						}
					}
					if ($names && $synonyms){
						$transcript->{attributes}->{'synonym'} = [@$names,@$synonyms];
					}
					elsif ($names){
						$transcript->{attributes}->{'synonym'} = $names;
					}
				}
				if ($tr_properties->{$tr_stable_id}{'synonym'}{'index'} && $tr_properties->{$tr_stable_id}{'synonym'}{'index'} <= 1){
					$transcript->{attributes}->{'display_name'} = $transcript->{attributes}->{'synonym'}->[0];
				}


				$transcript->{attributes}->{'stable_id'} = $tr_stable_id;
				$transcript->{attributes}->{'translation_stable_id'} = $tl_stable_id || $transcript->{attributes}->{'stable_id'};

			}
			if ($exception){
				push @exceptions,$gene;
			}
			else {
				push @valid,$gene;
			}

		}
		while (my $gene = shift @valid){
			if (my $out = $gene->structured_output(1)){
				print OUT $out;
			}
		}
		while (my $gene = shift @exceptions){
			if (my $out = $gene->structured_output(1)){
				print EXC $out;
			}
		}
	}

	if (-z $filename.".exception.gff.tmp"){
    unlink $filename.".exception.gff.tmp"
  }
  else {
    system "mv ".$filename.".exception.gff.tmp ".$filename.".exception.gff";
  }
	if (%dups){
		warn "WARNING: duplicate stable_ids introduced in [FILES] GFF\n";
		system "mv ".$filename.".gff.tmp ".$filename.".dup.gff";
		dedup_gff($filename,".dup.gff",\%dups);
	}
	system "mv ".$filename.".gff.tmp ".$filename.".gff" if -s $filename.".gff.tmp";
	return $filename.".gff" if -e $filename.".gff";

}

sub dedup_gff {
	my ($filename,$suffix,$dups) = @_;
	return if -z $filename.$suffix;
	unshift @ARGV,$filename.$suffix;
	open OUT, ">".$filename.".gff.tmp";
	my $gff = GFFTree->new({});
	$gff->name('GFF');
	my %counts;
	while ($gff->parse_chunk('change','region')){
		my @genes = $gff->by_type('gene');
		while (my $gene = shift @genes){
			if ($dups->{'gene'}{$gene->{attributes}->{'stable_id'}}){
				$counts{'gene'}{$gene->{attributes}->{'stable_id'}}++;
				$gene->{attributes}->{'stable_id'} .= '_'.$counts{'gene'}{$gene->{attributes}->{'stable_id'}};
			}
			my @daughters = $gene->daughters();
			while (my $transcript = shift @daughters){
				if ($dups->{'transcript'}{$transcript->{attributes}->{'stable_id'}}){
					$counts{'transcript'}{$transcript->{attributes}->{'stable_id'}}++;
					$transcript->{attributes}->{'stable_id'} .= '_'.$counts{'transcript'}{$transcript->{attributes}->{'stable_id'}};
				}
				if ($dups->{'translation'}{$transcript->{attributes}->{'translation_stable_id'}}){
					$counts{'translation'}{$transcript->{attributes}->{'translation_stable_id'}}++;
					$transcript->{attributes}->{'translation_stable_id'} .= '_'.$counts{'translation'}{$transcript->{attributes}->{'translation_stable_id'}};
				}
			}
			if (my $out = $gene->structured_output(1)){
				print OUT $out;
			}
		}
	}
	close OUT;
	if (%counts){
		foreach my $table (qw(gene transcript translation)){
			foreach my $stable_id (keys %{$counts{$table}}){
				warn "WARNING: $table stable_id=$stable_id count=".($counts{$table}{$stable_id}+1)."\n";
			}
		}
	}
	return;
}



sub get_properties {
	my ($infile,$property,$properties,$stable_ids_ref,$prop_ref) = @_;
	my ($index,$prop_location,$prop_regex,$prop_substitution) = @$prop_ref;
	my ($stable_id_location,$stable_id_regex,$stable_id_substitution) = @$stable_ids_ref;
	if ($infile->{'type'} eq 'fas'){
		# fasta file, read with Bio::Perl;
		my $seqin = Bio::SeqIO->new(-file => $infile->{'name'}, -format => 'fasta') || die "ERROR: unable to parse ".$infile->{'name'}." as fasta\n";
		while (my $seq = $seqin->next_seq()){
			my ($stable_id,$property_value);
			if ($stable_id_location eq 'DISPLAY_ID'){
				$stable_id = $seq->display_id();
			}
			elsif ($stable_id_location eq 'DESCRIPTION'){
				$stable_id = $seq->desc();
			}
			else {
				die "ERROR: could not recognise $stable_id_location as a location in a file of type ".$infile->{'type'}."\n";
			}
			if ($prop_location eq 'DISPLAY_ID'){
				$property_value = $seq->display_id();
			}
			elsif ($prop_location eq 'DESCRIPTION'){
				$property_value = $seq->desc();
			}
			else {
				die "ERROR: could not recognise $prop_location as a location in a file of type ".$infile->{'type'}."\n";
			}
			_match_property_to_stable_id($property,$properties,$index,$stable_id,$property_value,$stable_id_regex,$stable_id_substitution,$prop_regex,$prop_substitution);
		}
	}
	if ($infile->{'type'} =~ m/(csv|tsv)/){
		my $sep = $1 =~ m/csv/ ? ',' : '\t';
		open IN,$infile->{'name'} || die "ERROR: unable to open ".$infile->{'name'}."\n";
		$stable_id_location =~ s/FIELD_//;
		if ($stable_id_location !~ m/^\d+$/){
			die "ERROR: location in files of type ".$infile->{'type'}." should be formatted like 'FIELD_n', where n is the column number to select\n";
		}
		$prop_location =~ s/FIELD_//;
		if ($prop_location !~ m/^\d+$/){
			die "ERROR: location in files of type ".$infile->{'type'}." should be formatted like 'FIELD_n', where n is the column number to select\n";
		}
		while (<IN>){
			chomp;
			my @row = split /$sep/;
			my ($stable_id,$property_value);
			$stable_id = $row[($stable_id_location-1)] || die "ERROR FIELD_$stable_id_location is out of bounds in file ".$infile->{'name'}."\n";
			$property_value = $row[($prop_location-1)] || die "ERROR FIELD_$prop_location is out of bounds in file ".$infile->{'name'}."\n";
			_match_property_to_stable_id($property,$properties,$index,$stable_id,$property_value,$stable_id_regex,$stable_id_substitution,$prop_regex,$prop_substitution);
		}
	}
	return 1;
}


sub _match_property_to_stable_id {
	my ($property,$properties,$index,$stable_id,$property_value,$stable_id_regex,$stable_id_substitution,$property_regex,$property_substitution) = @_;
	$stable_id = _match_and_replace($stable_id,$stable_id_regex,$stable_id_substitution) || $stable_id;
	if ($property_value){
		$property_value = _match_and_replace($property_value,$property_regex,$property_substitution);
		if (!$properties->{$stable_id} || !$properties->{$stable_id}{$property} || $index < $properties->{$stable_id}{$property}{'index'}){
			# only update the description if there is not already a higher ranked description
			# would be more efficient to do this check earlier
			unshift @{$properties->{$stable_id}{$property}{'values'}},$property_value;
			$properties->{$stable_id}{$property}{'index'} = $index;
		}
		elsif ($properties->{$stable_id} && $properties->{$stable_id}{$property}){
			push @{$properties->{$stable_id}{$property}{'values'}},$property_value;
			$property_value = $properties->{$stable_id}{$property}{'values'}[0];
		}
	}
	elsif ($properties->{$stable_id} && $properties->{$stable_id}{$property}){
		$property_value = $properties->{$stable_id}{$property}{'values'}[0];
	}
	# TODO: check that values are not repeated in the array;
	if ($properties->{$stable_id}{$property}{'values'}){
		my @arr = @{$properties->{$stable_id}{$property}{'values'}};
		@{$properties->{$stable_id}{$property}{'values'}} = ();
		my %vals;
		for (my $i = 0; $i < @arr; $i++){
			if ($arr[$i] && !$vals{$arr[$i]}){
				push @{$properties->{$stable_id}{$property}{'values'}}, $arr[$i];
				$vals{$arr[$i]} = 1;
			}
		}


	}
	return ($stable_id,$property_value);
}


sub count_sequences {
	my ($params,$infiles,$dbh) = @_;
	if ($params->{'FILES'}{'SCAFFOLD'}){
		my $sth = $dbh->prepare("SELECT count(*)
									FROM seq_region
										JOIN coord_system
											ON seq_region.coord_system_id = coord_system.coord_system_id
									WHERE coord_system.name LIKE ".$dbh->quote('scaffold'));
		$sth->execute;
		if ($sth->rows > 0){
			my $count = $sth->fetchrow_arrayref()->[0];
			print "$count scaffold sequences in database\n";
		}
	}
	if ($params->{'FILES'}{'CONTIG'}){
		my $sth = $dbh->prepare("SELECT count(*)
									FROM seq_region
										JOIN coord_system
											ON seq_region.coord_system_id = coord_system.coord_system_id
									WHERE coord_system.name LIKE ".$dbh->quote('contig'));
		$sth->execute;
		if ($sth->rows > 0){
			my $count = $sth->fetchrow_arrayref()->[0];
			print "$count contig sequences in database\n";
		}
	}
	return 1;
}


sub add_seq_region_synonyms {
	my ($params,$infiles,$dbh) = @_;
	my $assembly_name = $params->{'META'}{'ASSEMBLY.NAME'};
	my @header;
	if ($infiles->{'SCAFFOLD_NAMES'}){
		open IN,$infiles->{'SCAFFOLD_NAMES'}{'name'};
		if ($params->{'SCAFFOLD_NAMES'}{'HEADER'}){
			<IN>;
			chomp;
			@header = split /\t/;
		}
		while (<IN>){
			chomp $_;
			my @names = split /\t/,$_;
			shift @names if ($names[0] =~ m/^$assembly_name/);
			unless (_add_seq_synonyms($dbh,undef,@names)){
     		    warn "WARNING: could not find any matching seq_region.name or seq_region_synonym.synonym on line $.\n";
			}
		}
	}
	if ($params->{'SCAFFOLD_NAMES'}{'SCAFFOLD'}){
		my $sth = $dbh->prepare("SELECT seq_region.seq_region_id, seq_region.name
									FROM seq_region
										JOIN coord_system
											ON seq_region.coord_system_id = coord_system.coord_system_id
									WHERE coord_system.name LIKE 'scaffold'");
		$sth->execute;
		while (my($seq_region_id,$seq_region_name) = $sth->fetchrow_array()){
			my $name = _match_and_replace($seq_region_name,$params->{'SCAFFOLD_NAMES'}{'SCAFFOLD'}[0],$params->{'SCAFFOLD_NAMES'}{'SCAFFOLD'}[1]);
			_add_seq_synonyms($dbh,$seq_region_id,$name);
		}
	}
	if ($params->{'SCAFFOLD_NAMES'}{'CONTIG'}){
		my $sth = $dbh->prepare("SELECT seq_region.seq_region_id, seq_region.name
									FROM seq_region
										JOIN coord_system
											ON seq_region.coord_system_id = coord_system.coord_system_id
									WHERE coord_system.name LIKE 'contig'");
		$sth->execute;
		while (my($seq_region_id,$seq_region_name) = $sth->fetchrow_array()){
			my $name = _match_and_replace($seq_region_name,$params->{'SCAFFOLD_NAMES'}{'CONTIG'}[0],$params->{'SCAFFOLD_NAMES'}{'CONTIG'}[1]);
			_add_seq_synonyms($dbh,$seq_region_id,$name);
		}
	}
	return 1;
}


sub _match_and_replace {
	my ($str,$match,$replace) = @_;
	if ($match){
		$match =~ s/^\///;
		$match =~ s/\/$//;
		if (ref $str ne 'ARRAY'){
			$str =~ m/$match/;
			$str = $1;
		}
		else {
			foreach my $element (@$str){
				if ($element =~ m/$match/){
					$str = $1;
					last;
				}
			}
		}
		if ($str && $replace){
			$replace =~ s/^\///;
			$replace =~ s/\/$//;
			$replace =~ s/\\\//__,,__/;
			my ($sub,$with) = split /\//,$replace;
			$sub =~ s/__,,__/\//;
			$with =~ s/__,,__/\//;
			if (!$with){
				$str =~ s/$sub//;
			}
			elsif ($with !~ /\$/){
				$str =~ s/$sub/$with/;
			}
			else {
				my ($left,undef,$right) = split /(\$1)/,$with;
				$left ||= '';
				$right ||= '';
				if ($str =~ m/$sub/){
					my $new_value = $left.$1.$right;
					$str =~ s/$sub/$new_value/;
				}
			}
		}
	}
	else {
		$str = $str->[0];
	}
	return $str;
}

sub _add_seq_synonyms {
	my $dbh = shift;
	my $seq_region_id = shift;
	if (!$seq_region_id){
		my $condition = "name LIKE '".join("' OR name LIKE '",@_)."'";
		$condition =~ s/name\sLIKE\s'na'\sOR\s//;
		$condition =~ s/OR\name\sLIKE\s'na'//;
		my $sth1 = $dbh->prepare("SELECT seq_region_id FROM seq_region WHERE $condition");
		$sth1->execute;
		if ($sth1->rows > 0){
			$seq_region_id = $sth1->fetchrow_arrayref()->[0];
		}
		unless ($seq_region_id){
			$condition =~ s/name/synonym/g;
			my $sth2 = $dbh->prepare("SELECT seq_region_id FROM seq_region_synonym WHERE $condition");
    		$sth2->execute;
			if ($sth2->rows > 0){
				$seq_region_id = $sth2->fetchrow_arrayref()->[0];
			}
		}
		return unless $seq_region_id;
	}
	while (my $region = shift @_){
		next if $region eq 'na';
		my $sth3 = $dbh->prepare("SELECT seq_region_id FROM seq_region WHERE seq_region_id = $seq_region_id AND name LIKE ".$dbh->quote($region));
		$sth3->execute;
		unless ($sth3->rows > 0){
			my $sth4 = $dbh->prepare("SELECT seq_region_synonym_id FROM seq_region_synonym WHERE seq_region_id = $seq_region_id AND synonym LIKE ".$dbh->quote($region));
			$sth4->execute;
			unless ($sth4->rows > 0){
				$dbh->do("INSERT INTO seq_region_synonym (seq_region_id, synonym) VALUES ($seq_region_id, ".$dbh->quote($region).")");
			}
		}
	}
	return 1;
}


sub load_sequences {
	warn "trying";
	my ($params,$infiles,$dbh) = @_;
	my $connection_info = 	 '-dbuser '.$params->{'DATABASE_CORE'}{'RW_USER'}
							.' -dbhost '.$params->{'DATABASE_CORE'}{'HOST'}
							.' -dbport '.$params->{'DATABASE_CORE'}{'PORT'}
							.' -dbname '.$params->{'DATABASE_CORE'}{'NAME'}
							.' -dbpass '.$params->{'DATABASE_CORE'}{'RW_PASS'};
	my $perl_libs = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl/modules';
	my $perl = "perl -I $perl_libs";
	if ($infiles->{'CONTIG'} && !$infiles->{'SCAFFOLD'}){
		if ($infiles->{'CONTIG'}{'type'} ne 'fas'){
			die "ERROR: [FILES] CONTIG $infiles->{'CONTIG'}{'name'} must be a fasta file if it is the only sequence file specified\n";
		}
		system $perl.' '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline/scripts/load_seq_region.pl '.$connection_info.' -coord_system_name contig -rank 1 -coord_system_version '.$params->{'META'}{'ASSEMBLY.NAME'}.' -default_version -sequence_level -verbose -fasta_file '.$infiles->{'CONTIG'}{'name'}.' -replace_ambiguous_bases';
	}
	elsif ($infiles->{'CONTIG'}){
		if ($infiles->{'CONTIG'}{'type'} ne 'fas' && $infiles->{'SCAFFOLD'}{'type'} ne 'fas'){
			die "ERROR: at least one of [FILES] CONTIG $infiles->{'CONTIG'}{'name'} or SCAFFOLD $infiles->{'SCAFFOLD'}{'name'} must be a fasta file\n";
		}
		if ($infiles->{'CONTIG'}{'type'} ne 'agp' && $infiles->{'SCAFFOLD'}{'type'} ne 'agp'){
			die "ERROR: at least one of [FILES] CONTIG $infiles->{'CONTIG'}{'name'} or SCAFFOLD $infiles->{'SCAFFOLD'}{'name'} must be agp file\n";
		}
		if ($infiles->{'CONTIG'}{'type'} ne 'agp'){
			system $perl.' '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline/scripts/load_seq_region.pl '.$connection_info.' -coord_system_name scaffold -rank 1 -coord_system_version '.$params->{'META'}{'ASSEMBLY.NAME'}.' -default_version -agp_file '.$infiles->{'SCAFFOLD'}{'name'};
			system $perl.' '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline/scripts/load_seq_region.pl '.$connection_info.' -coord_system_name contig -rank 2 -default_version -sequence_level -verbose -fasta_file '.$infiles->{'CONTIG'}{'name'}.' -replace_ambiguous_bases';
			system $perl.' '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline/scripts/load_agp.pl '.$connection_info.' -assembled_name scaffold -assembled_version '.$params->{'META'}{'ASSEMBLY.NAME'}.' -component_name contig -agp_file '.$infiles->{'SCAFFOLD'}{'name'};
		}
		else {
			system $perl.' '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline/scripts/load_seq_region.pl '.$connection_info.' -coord_system_name scaffold -rank 1 -default_version -sequence_level -verbose -fasta_file '.$infiles->{'SCAFFOLD'}{'name'}.' -replace_ambiguous_bases';
			system $perl.' '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline/scripts/load_seq_region.pl '.$connection_info.' -coord_system_name contig -rank 2 -coord_system_version '.$params->{'META'}{'ASSEMBLY.NAME'}.' -default_version -agp_file '.$infiles->{'CONTIG'}{'name'};
			system $perl.' '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline/scripts/load_agp.pl '.$connection_info.' -assembled_name contig -assembled_version '.$params->{'META'}{'ASSEMBLY.NAME'}.' -component_name scaffold -agp_file '.$infiles->{'CONTIG'}{'name'};
		}
		$dbh->{RaiseError} = 0;
		$dbh->do('INSERT INTO meta(species_id, meta_key,meta_value) VALUES (1, '.$dbh->quote('assembly.mapping').','.$dbh->quote('scaffold:'.$params->{'META'}{'ASSEMBLY.NAME'}.'|contig').')');
		$dbh->{RaiseError} = 1;
	}
	else {
		if ($infiles->{'SCAFFOLD'}{'type'} ne 'fas'){
			die "ERROR: [FILES] SCAFFOLD $infiles->{'SCAFFOLD'}{'name'} must be a fasta file if it is the only sequence file specified\n";
		}
		system $perl.' '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline/scripts/load_seq_region.pl '.$connection_info.' -coord_system_name scaffold -rank 1 -coord_system_version '.$params->{'META'}{'ASSEMBLY.NAME'}.' -default_version -sequence_level -verbose -fasta_file '.$infiles->{'SCAFFOLD'}{'name'}.' -replace_ambiguous_bases';

	}
	system $perl.' '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline/scripts/set_toplevel.pl '.$connection_info;

	return 1;

}

sub setup_core_db {
	my $dbh = shift;
	my $params = shift;

	$dbh->do('DROP DATABASE IF EXISTS '.$params->{'DATABASE_CORE'}{'NAME'}) || die "ERROR: unable to drop existing [DATABASE_CORE] NAME ".$params->{'DATABASE_CORE'}{'NAME'}." using provided settings";
	$dbh->do('CREATE DATABASE '.$params->{'DATABASE_CORE'}{'NAME'}) || die "ERROR: unable to create to [DATABASE_CORE] NAME ".$params->{'DATABASE_CORE'}{'NAME'}." using provided settings";

	my $file = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl/sql/table.sql';
	my $connection_info = 	 '-u'.$params->{'DATABASE_CORE'}{'RW_USER'}
							.' -h'.$params->{'DATABASE_CORE'}{'HOST'}
							.' -P'.$params->{'DATABASE_CORE'}{'PORT'}
							.' -D'.$params->{'DATABASE_CORE'}{'NAME'}
							.' -p'.$params->{'DATABASE_CORE'}{'RW_PASS'};
	system "mysql $connection_info < $file";

	$file = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline/sql/table.sql';
	system "mysql $connection_info < $file";

	$dbh->do('USE '.$params->{'DATABASE_CORE'}{'NAME'});

	## loop through [META] and load data
	$params->{'META'}{'GENEBUILD.LEVEL'} ||= 'toplevel';
	$params->{'META'}{'TRANSCRIPTBUILD.LEVEL'} ||= 'toplevel';
	$params->{'META'}{'EXONBUILD.LEVEL'} ||= 'toplevel';
	foreach my $key (sort keys %{$params->{'META'}}){
		my $value = $params->{'META'}{$key};
		$key = lc $key;
		$key = 'species.'.$key unless $key =~ m/\./;
		if (ref $value){
			for (my $i = 0; $i < @{$value}; $i++){
				my $tmp = $value->[$i];
				$dbh->do('INSERT INTO meta (meta_key,meta_value) VALUES ('.$dbh->quote($key).','.$dbh->quote($tmp).')');
			}
		}
		else {
			$dbh->do('INSERT INTO meta (meta_key,meta_value) VALUES ('.$dbh->quote($key).','.$dbh->quote($value).')');
		}
	}


	## add meta_key = taxonomy
	my $dbh_tax = taxonomy_db_connect($params);
	my $sth_tax = $dbh_tax->prepare("SELECT names.taxonid, names.name FROM names JOIN nodes ON names.taxonid = nodes.parenttaxonid WHERE nodes.taxonid = ? limit 1");
	$sth_tax->execute($params->{'META'}{'SPECIES.TAXONOMY_ID'});
	my $last_taxid = -1;
	while (my ($taxid,$name) = $sth_tax->fetchrow_array()){
		last if $taxid == $last_taxid;
		$dbh->do("INSERT INTO meta (meta_key, meta_value) VALUES ('species.classification', '$name')");
  		$sth_tax->execute($taxid);
  		$last_taxid = $taxid;
  	}

	## add interpro table
	system 'mysqldump --single-transaction'
			.' --user='.$params->{'DATABASE_TEMPLATE'}{'RO_USER'}
			.' --host='.$params->{'DATABASE_TEMPLATE'}{'HOST'}
			.' --port='.$params->{'DATABASE_TEMPLATE'}{'PORT'}
			.' '.$params->{'DATABASE_TEMPLATE'}{'NAME'}
			.' interpro'
			.' | mysql '
			.' -h '.$params->{'DATABASE_CORE'}{'HOST'}
			.' -P '.$params->{'DATABASE_CORE'}{'PORT'}
			.' -u '.$params->{'DATABASE_CORE'}{'RW_USER'}
			.' -p'.$params->{'DATABASE_CORE'}{'RW_PASS'}
			.' '.$params->{'DATABASE_CORE'}{'NAME'};



	## populate production db tables
	system 	 'perl '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-production/scripts/production_database/populate_production_db_tables.pl'
			.' -h '.$params->{'DATABASE_CORE'}{'HOST'}
			.' -P '.$params->{'DATABASE_CORE'}{'PORT'}
			.' -u '.$params->{'DATABASE_CORE'}{'RW_USER'}
			.' -p '.$params->{'DATABASE_CORE'}{'RW_PASS'}
			.' -d '.$params->{'DATABASE_CORE'}{'NAME'}
			.' -mh '.$params->{'DATABASE_TEMPLATE'}{'HOST'}
			.' -mP '.$params->{'DATABASE_TEMPLATE'}{'PORT'}
			.' -mu '.$params->{'DATABASE_TEMPLATE'}{'RO_USER'}
			.' -md '.$params->{'DATABASE_TEMPLATE'}{'NAME'}
			.' -dp ./';

	## clean up sql backups
	system 'rm '.$params->{'DATABASE_CORE'}{'NAME'}.'*.sql';

	return 1;
}

sub agp_file_summary {
	my ($params,$infile,$type) = @_;
	my ($scaffold_hash_ref,%all_stats,%N50nums);
	$scaffold_hash_ref = _agp_file_to_hash($infile->{'name'},$type);
	$all_stats{$type} = _seq_hash_to_stats($scaffold_hash_ref->{$type},'len');
    $N50nums{$type} = $all_stats{$type}{num_n50};
    $all_stats{'contig'} = _seq_hash_to_stats($scaffold_hash_ref->{'contig'},'len');
    $all_stats{'Ns'} = _seq_hash_to_stats($scaffold_hash_ref->{'n_hash'},'len');
    mkdir "summary";
	open OUT,">summary/$type.stats.txt";
    print OUT ucfirst(lc($type)).":\n";
    print OUT "LongestScaffold\t".$all_stats{$type}{max_length},"\n";
	print OUT "Num\t".$all_stats{$type}{num_seqs},"\n";
    print OUT "Span\t".$all_stats{$type}{span},"\n";
    print OUT "Min\t".$all_stats{$type}{min_length},"\n";
    print OUT "Mean\t".$all_stats{$type}{mean_length},"\n";
    print OUT "N50\t".$all_stats{$type}{n50},"\n";
    print OUT "NumN50\t".$all_stats{$type}{num_n50},"\n";
    if ($all_stats{'contig'}{num_seqs} && $all_stats{'contig'}{num_seqs} > 0){
    	$N50nums{'CONTIG'} = $all_stats{'contig'}{num_n50};
    	print OUT "Contigs:\n";
    	print OUT "Num\t".$all_stats{'contig'}{num_seqs},"\n";
    	print OUT "Span\t".$all_stats{'contig'}{span},"\n";
    	print OUT "N50\t".$all_stats{'contig'}{n50},"\n";
    	print OUT "NumN50\t".$all_stats{'contig'}{num_n50},"\n";
    }
	if ($all_stats{'Ns'}{num_seqs} && $all_stats{'Ns'}{num_seqs} > 0){
    	$N50nums{'Ns'} = $all_stats{'Ns'}{num_n50};
    	print OUT "Runs of >=10 Ns:\n";
    	print OUT "Num\t".$all_stats{'Ns'}{num_seqs},"\n";
    	print OUT "Span\t".$all_stats{'Ns'}{span},"\n";
    	print OUT "N50\t".$all_stats{'Ns'}{n50},"\n";
    	print OUT "NumN50\t".$all_stats{'Ns'}{num_n50},"\n";
    }
    close OUT;
    open RDATA,">summary/$type.stats.Rdata";
    map { print RDATA "Scaffold\t$_\n" } @{$all_stats{$type}{lengths}};
    map { print RDATA "Contig\t$_\n" }   @{$all_stats{'contig'}{lengths}};
    map { print RDATA "N\t$_\n" }        @{$all_stats{'Ns'}{lengths}};
    close RDATA;
	make_cumulative_plot("summary/$type.stats.Rdata",\%N50nums);
	my %return_stats;
	$return_stats{'assembly_span'} = $all_stats{$type}{span};
	$return_stats{'lengths'} = $all_stats{'SCAFFOLD'}{lengths} if $all_stats{'SCAFFOLD'};
	$return_stats{'contig_lengths'} = $all_stats{'contig'}{lengths} if $all_stats{'contig'};
	$return_stats{'assembly_atgc'} = $return_stats{'assembly_span'};
	$return_stats{'assembly_atgc'} -= $all_stats{'Ns'}{span} if ($all_stats{'Ns'}{num_seqs} && $all_stats{'Ns'}{num_seqs} > 0);
	$return_stats{lc($type).'_count'} = $all_stats{$type}{num_seqs};
	return \%return_stats;
}

sub _agp_file_to_hash {
    my ($agpfile,$type) = @_;
    my (%main_hash,%scaf_hash,%contig_hash,%n_hash);
    open IN,$agpfile;
    my $n_id = 0;
    while (<IN>){
    	next if m/^#/;
    	chomp;
    	my ($seqid,$start,$end,$ctr,$flag,$contig) = split /\t/;
    	$scaf_hash{$seqid}{seq} += $end - $start + 1;
    	if ($flag eq 'W'){
    		$contig_hash{$contig}{seq} += $end - $start + 1
        }
        elsif ($flag eq 'N'){
    		$n_hash{++$n_id}{seq} += $end - $start + 1
        }
    }
    $main_hash{$type} = \%scaf_hash;
    $main_hash{'contig'} = \%contig_hash;
    $main_hash{'n_hash'} = \%n_hash;
    return \%main_hash;
}

sub fasta_file_summary {
	my ($params,$infile,$type) = @_;
	my (%scaffold_hash,%all_stats,%N50nums);
	$scaffold_hash{$type} = _fasta_file_to_scaffold_hash($infile->{'name'});
	$scaffold_hash{'contig'} = _scaffold_hash_to_contig_hash($scaffold_hash{$type});
	$scaffold_hash{'n_hash'} = _scaffold_hash_to_n_hash($scaffold_hash{$type});
	$all_stats{$type} = _seq_hash_to_stats($scaffold_hash{$type},'seq');
    $N50nums{$type} = $all_stats{$type}{num_n50};
    $all_stats{'contig'} = _seq_hash_to_stats($scaffold_hash{'contig'},'seq');
    $all_stats{'Ns'} = _seq_hash_to_stats($scaffold_hash{'n_hash'},'seq');


    mkdir "summary";
	open OUT,">summary/$type.stats.txt";
    print OUT ucfirst(lc($type)).":\n";
    print OUT "LongestScaffold\t".$all_stats{$type}{max_length},"\n";
	print OUT "Num\t".$all_stats{$type}{num_seqs},"\n";
    print OUT "Span\t".$all_stats{$type}{span},"\n";
    print OUT "Min\t".$all_stats{$type}{min_length},"\n";
    print OUT "Mean\t".$all_stats{$type}{mean_length},"\n";
    print OUT "N50\t".$all_stats{$type}{n50},"\n";
    print OUT "NumN50\t".$all_stats{$type}{num_n50},"\n";
    print OUT "GC\t".$all_stats{$type}{gc},"\n";
    if ($type ne 'CONTIG' && $all_stats{'contig'}{num_seqs} && $all_stats{'contig'}{num_seqs} > 0){
    	$N50nums{'CONTIG'} = $all_stats{'contig'}{num_n50};
    	print OUT "Contigs:\n";
    	print OUT "Num\t".$all_stats{'contig'}{num_seqs},"\n";
    	print OUT "Span\t".$all_stats{'contig'}{span},"\n";
    	print OUT "N50\t".$all_stats{'contig'}{n50},"\n";
    	print OUT "NumN50\t".$all_stats{'contig'}{num_n50},"\n";
    }
	if ($all_stats{'Ns'}{num_seqs} && $all_stats{'Ns'}{num_seqs} > 0){
    	$N50nums{'Ns'} = $all_stats{'Ns'}{num_n50};
		print OUT "Runs of >=10 Ns:\n";
    	print OUT "Num\t".$all_stats{'Ns'}{num_seqs},"\n";
    	print OUT "Span\t".$all_stats{'Ns'}{span},"\n";
    	print OUT "N50\t".$all_stats{'Ns'}{n50},"\n";
    	print OUT "NumN50\t".$all_stats{'Ns'}{num_n50},"\n";
    }
    open RDATA,">summary/$type.stats.Rdata";
    map { print RDATA ucfirst(lc($type))."\t$_\n" } @{$all_stats{$type}{lengths}};
    if ($type ne 'CONTIG' && $all_stats{'contig'}{num_seqs} && $all_stats{'contig'}{num_seqs} > 0){
    	map { print RDATA "Contig\t$_\n" } @{$all_stats{'contig'}{lengths}};
    }
    if ($all_stats{'Ns'}{num_seqs} && $all_stats{'Ns'}{num_seqs} > 0){
    	map { print RDATA "N\t$_\n" } @{$all_stats{'Ns'}{lengths}};
    }
    close RDATA;
	make_cumulative_plot("summary/$type.stats.Rdata",\%N50nums);

	my %return_stats;
	$return_stats{'assembly_span'} = $all_stats{$type}{span};
	$return_stats{'assembly_atgc'} = $return_stats{'assembly_span'};
	$return_stats{'assembly_atgc'} -= $all_stats{'Ns'}{span} if ($all_stats{'Ns'}{num_seqs} && $all_stats{'Ns'}{num_seqs} > 0);
	$return_stats{'assembly_gc'} = $all_stats{$type}{gc};
	$return_stats{'lengths'} = $all_stats{'SCAFFOLD'}{lengths} if $all_stats{'SCAFFOLD'};
	$return_stats{'contig_lengths'} = $all_stats{'contig'}{lengths} if $all_stats{'contig'};
	$return_stats{'GCs'} = $all_stats{'SCAFFOLD'}{'GCs'} if $all_stats{'SCAFFOLD'}{'GCs'};
	$return_stats{'Ns'} = $all_stats{'SCAFFOLD'}{'Ns'} if $all_stats{'SCAFFOLD'}{'Ns'};
	$return_stats{lc($type).'_count'} = $all_stats{$type}{num_seqs};
	return \%return_stats;
}


sub _fasta_file_to_scaffold_hash {
    my $fastafile = shift;
    my %sequences;
    open IN,$fastafile;
    my $seqid;
    while (<IN>){
        next if /^\s*$/ or /^#/;
		if (/^>\s*(\S+)/) {
            $seqid = $1;
        }
        else {
            chomp($sequences{$seqid}{seq} .= $_  );
        }
    }
    return \%sequences;
}

sub _scaffold_hash_to_contig_hash {
	my $scaffoldhash = shift;
    my %contighash;
    my $contigid = 0;
    my $numNtosplitscaffold = 10;
    for my $seqid (keys %{$scaffoldhash}) {
        my $seq   = $$scaffoldhash{$seqid}{seq};
        foreach (split (/N{$numNtosplitscaffold,}/i,$seq)) {
            next if /^\s*$/; # skip blank entries which you get if seq begins with Ns
            next if length($_) < 100;
            $contighash{++$contigid}{seq} = $_;
        }
    }
    return \%contighash;
}

sub _scaffold_hash_to_n_hash {
	my $scaffoldhash = shift;
    my %n_hash;
    my $n_id = 0;
    for my $seqid (keys %{$scaffoldhash}) {
        my $seq = $$scaffoldhash{$seqid}{seq};
        foreach (split (/[atgc]+/i,$seq)) {
            next if /^\s*$/; # skip blank entries which you get if seq begins with Ns
            next if length($_) < 10;
            $n_hash{++$n_id}{seq} = $_;
        }
    }
    return \%n_hash;
}


sub _seq_hash_to_stats {
	my ($seqhash,$datatype) = @_;
    my %stats;
    my $gc_sum = 0;
    my $atgc_sum = 0;
    my @sorted_lengths;
    my @sorted_seqs;
    foreach (keys %{$seqhash} ) {
        my $seq = $$seqhash{$_}{seq};
        if ($datatype eq 'seq'){
	        push @sorted_lengths, length ($seq);
    	    push @sorted_seqs, $seq;
    	    $gc_sum               += ($seq =~ tr/gcGC/gcGC/);
        	$atgc_sum             += ($seq =~ tr/atgcATGC/atgcATGC/);
        }
        else {
        	push @sorted_lengths, $seq;
        }
    }
    $stats{gc}            = $atgc_sum        ? sprintf("%.3f",$gc_sum/$atgc_sum)   : 0;
    @sorted_lengths       = sort {$b <=> $a} @sorted_lengths;
    $stats{lengths}       = \@sorted_lengths;
    $stats{min_length}    =  @sorted_lengths ? $sorted_lengths[$#sorted_lengths]   : 0;
    $stats{max_length}    =  @sorted_lengths ? $sorted_lengths[0]                  : 0;
    $stats{num_seqs}      =  @sorted_lengths ? @sorted_lengths                     : 0;
    $stats{span}          =  @sorted_lengths ? sum(@sorted_lengths)                : 0;
    $stats{mean_length}   =  @sorted_lengths ? int($stats{span}/$stats{num_seqs})  : 0;
    $stats{median_length} =  @sorted_lengths ? $sorted_lengths[$stats{num_seqs}/2] : 0;
    my ($csum, $nlen, $n50, $num_n50) = (0,0,0,0);
    for $nlen (@sorted_lengths) { $csum += $nlen; $num_n50++; $n50 = $nlen; last if $csum >= ($stats{span}/2) }
    $stats{n50}           = $n50;
    $stats{num_n50}       = $num_n50;
    my $seq = '';
    my $thousandth = int($stats{span} / 1000);
    my @GCs;
	my @Ns;
	while (my $nextseq = shift @sorted_seqs) {
		$seq .= $nextseq;
		while (length $seq >= $thousandth){
			my $chunk = substr($seq,0,$thousandth);
			my $n = () = $chunk =~ /n/gi;
			my $gc = () = $chunk =~ /[gc]/gi;
			$gc = ($thousandth - $n) ? $gc / ($thousandth - $n) : 0;
			push @Ns, int(( $n / $thousandth )*1000)/10;
			push @GCs, int($gc*1000)/10;
			$seq = substr($seq,$thousandth)
		}
	}
	while (length $seq >= $thousandth){
		my $chunk = substr($seq,0,$thousandth);
		my $n = () = $chunk =~ /n/gi;
		my $gc = () = $chunk =~ /[gc]/gi;
		$gc = ($thousandth - $n) ? $gc / ($thousandth - $n) : 0;
		push @Ns, int(( $n / $thousandth )*1000)/10;
		push @GCs, int($gc*1000)/10;
		$seq = substr($seq,$thousandth)
	}
	$stats{GCs} = @GCs ? \@GCs : undef;
    $stats{Ns} = @Ns ? \@Ns : undef;
    return \%stats;
}

sub make_cumulative_plot {

my ($rdatafile,$N50nums) = @_;
open RSCRIPT, ">$rdatafile.R" or die $!;

print RSCRIPT '
args <- commandArgs(trailingOnly = TRUE)
data=read.delim(args[1],header=F,col.names=c("type","len"))
for (t in levels(data$type)) {
     data[data$type==t,"cumsum"] =cumsum   (data[data$type==t,"len"])
     data[data$type==t,"Rank"]   =seq_along(data[data$type==t,"len"])
}
cols = c(rgb(166/255,206/255,227/255),rgb(31/255,120/255,180/255),rgb(178/255,223/255,138/255))
xmax = max(data[,"Rank"])
ymax = max(data[,"cumsum"])
lwd=2
png(paste(args[1],"png",sep="."),900,900,res=200)
par(las=1,mar=c(4,5,1,1),cex.axis=0.8,cex.lab=1)
plot(1,1,type="n",xlab="",ylab="",xlim=c(0,xmax),ylim=c(0,ymax))
grid(col = "lightgray",lty=1, lwd = 0.5, equilogs = TRUE)
mtext(text="Cumulative length (bp)", side=2, line=3.5, las=0, cex=1)
mtext(text="Sequences ranked by size", side=1, line=2.5, las=0, cex=1)
labels = c(NA)
legendcols = c(NA)
';
if ($N50nums->{'Ns'}){
print RSCRIPT 'n50 <- data[which(data$Rank=='.$N50nums->{'Ns'}.' & data$type=="N"),]
lines(c(n50$Rank,-xmax),c(n50$cumsum,n50$cumsum),col=cols[3],lwd=1,lty=2)
lines(c(n50$Rank,n50$Rank),c(n50$cumsum,-ymax),col=cols[3],lwd=1,lty=2)
labels = c("Ns",labels)
legendcols = c(cols[3],legendcols)
';
}
if ($N50nums->{'CONTIG'}){
print RSCRIPT 'n50 <- data[which(data$Rank=='.$N50nums->{'CONTIG'}.' & data$type=="Contig"),]
lines(c(n50$Rank,-xmax),c(n50$cumsum,n50$cumsum),col=cols[2],lwd=1,lty=2)
lines(c(n50$Rank,n50$Rank),c(n50$cumsum,-ymax),col=cols[2],lwd=1,lty=2)
labels = c("Contigs",labels)
legendcols = c(cols[2],legendcols)
';
}
if ($N50nums->{'SCAFFOLD'}){
print RSCRIPT 'n50 <- data[which(data$Rank=='.$N50nums->{'SCAFFOLD'}.' & data$type=="Scaffold"),]
lines(c(n50$Rank,-xmax),c(n50$cumsum,n50$cumsum),col=cols[1],lwd=1,lty=2)
lines(c(n50$Rank,n50$Rank),c(n50$cumsum,-ymax),col=cols[1],lwd=1,lty=2)
labels = c("Scaffolds",labels)
legendcols = c(cols[1],legendcols)
';
}
if ($N50nums->{'Ns'}){
print RSCRIPT 'lines(data[data$type=="N",c("Rank","cumsum")],type="l",lwd=lwd,col=cols[3])
';
}
if ($N50nums->{'CONTIG'}){
print RSCRIPT 'lines(data[data$type=="Contig",c("Rank","cumsum")],type="l",lwd=lwd,col=cols[2])
';
}
if ($N50nums->{'SCAFFOLD'}){
print RSCRIPT 'lines(data[data$type=="Scaffold",c("Rank","cumsum")],type="l",lwd=lwd,col=cols[1])
';
}
print RSCRIPT 'labels <- labels[-length(labels)]
legend("bottomright",labels,col=legendcols,lty=1,lwd=lwd,ncol=1,box.col = "black",bg = "white")
dev.off()
';
close RSCRIPT;
system "Rscript $rdatafile.R $rdatafile >/dev/null";
#unlink  "$rdatafile.R";
}

sub about_page {
	my $params = shift;
	my $production_name = $params->{'META'}{'SPECIES.PRODUCTION_NAME'};
	my $species = $params->{'META'}{'SPECIES.SCIENTIFIC_NAME'};
	open ABOUT,">summary/about_$production_name.html";
print ABOUT <<EOF;
<!-- {about} --><a name="about"></a>
  <h2 id="about" class="first">About <em>$species</em></h2>
  <p id="about"></p>
<!-- {about} -->
<!-- {assembly} --><a name="assembly"></a>
  <h2 id="assembly" style="margin-bottom:5px">Assembly</h2>
  <img style="width:100%;" src="/i/species/about/${production_name}_feature.length.png" title="distribution of feature lengths in the assembly"/>
  <br/>
  <img style="width:100%;" src="/i/species/about/${production_name}_feature.dist.png" title="distribution of features across scaffolds"/>
<!-- {assembly} -->
<!-- {annotation} --><a name="annotation"></a>
  <h2 id="annotation" style="margin-bottom:5px">Annotation</h2>
  <p></p>
<!-- {annotation} -->
<!-- {references} -->
  <h2>References</h3>
  <ol>
    <!--li><a id="ref-01"></a><a href="http://europepmc.org/abstract/MED/123456">TITLE</a>.<br>AUTHOR A, AUTHOR B. YYYY. JOURNAL. VOLUME:PAGES.</li-->
  </ol>
 <!-- {references} -->
EOF

close ABOUT;
	return 1;
}

sub stats_page {
	my ($params,$stats) = @_;;
	my $production_name = $params->{'META'}{'SPECIES.PRODUCTION_NAME'};
	my $species = $params->{'META'}{'SPECIES.SCIENTIFIC_NAME'};
	my $assembly_name = $params->{'META'}{'ASSEMBLY.NAME'};
	my $provider_name = $params->{'META'}{'PROVIDER.NAME'};
	my $provider_url = $params->{'META'}{'PROVIDER.URL'};
	my $assembly_date = $params->{'META'}{'ASSEMBLY.DATE'};
	my $db_version = $params->{'ENSEMBL'}{'BRANCH'}.".1";
	$db_version =~ s/release\///;
	my $assembly_span = $stats->{'assembly_span'};
	$assembly_span = readable($assembly_span,'bp');
	my $assembly_atgc = $stats->{'assembly_atgc'};
	$assembly_atgc = readable($assembly_atgc,'bp');
	my $assembly_gc_pct = $stats->{'assembly_gc'} * 100;
	my $genebuild_version = $params->{'META'}{'GENEBUILD.VERSION'};
	my $genebuild_method = $params->{'META'}{'GENEBUILD.METHOD'};
	my $gene_count = $stats->{'gene_count'};
	my $transcript_count = $stats->{'transcript_count'};
	my $translation_count = $stats->{'translation_count'};
	my $scaffold_count = $stats->{'scaffold_count'} || undef;
	my $contig_count = $stats->{'contig_count'} || undef;
	open STATS,">summary/stats_$production_name.html";
print STATS <<EOF;
<h3>Summary</h3>
<table class="ss tint species-stats">
  <tbody>
    <tr class="bg2">
      <td class="data">Assembly:</td>
      <td class="value">$assembly_name, <a href="$provider_url">$provider_name</a>, $assembly_date</td>
    </tr>
    <tr>
      <td class="data">Database version:</td>
      <td class="value">$db_version</td>
    </tr>
    <tr class="bg2">
      <td class="data">Base pairs:</td>
      <td class="value">$assembly_atgc</td>
    </tr>
    <tr>
      <td class="data">Golden path length:</td>
      <td class="value">$assembly_span</td>
    </tr>
    <tr class="bg2">
      <td class="data">GC content:</td>
      <td class="value">$assembly_gc_pct %</td>
    </tr>
    <tr>
      <td class="data">Genebuild version:</td>
      <td class="value">$genebuild_version</td>
    </tr>
    <tr class="bg2">
      <td class="data">Genebuild method:</td>
      <td class="value">$genebuild_method</td>
    </tr>
  </tbody>
</table>
<a href="/$production_name/Info/IPtop500?db=core">Table of top 500 InterPro hits</a><br/><br/>

<h3>Gene counts</h3>
<table class="ss tint species-stats">
  <tbody>
    <tr class="bg2">
      <td class="data"><span class="glossary_mouseover">Genes</span>:</td>
      <td class="value">$gene_count</td>
    </tr>
    <tr>
      <td class="data"><span class="glossary_mouseover">Transcripts</span>:</td>
      <td class="value">$transcript_count</td>
    </tr>
    <tr class="bg2">
      <td class="data"><span class="glossary_mouseover">Translations</span>:</td>
      <td class="value">$translation_count</td>
    </tr>
  </tbody>
</table>

<h3>Coordinate systems</h3>
<table class="ss tint species-stats">
  <tbody>
EOF
if ($scaffold_count){
print STATS <<EOF;
    <tr class="bg2">
      <td class="data">scaffold</td>
      <td class="value">$scaffold_count sequences</td>
    </tr>
EOF
}
if ($contig_count){
print STATS <<EOF;
    <tr>
      <td class="data">contig</td>
      <td class="value">$contig_count sequences</td>
    </tr>
EOF
}
print STATS <<EOF;
  </tbody>
</table>
EOF

	close STATS;
	return 1;
}


sub assembly_page {
	my ($params,$stats) = @_;
	my $production_name = $params->{'META'}{'SPECIES.PRODUCTION_NAME'};
	my $species = $params->{'META'}{'SPECIES.SCIENTIFIC_NAME'};
	my $span = $stats->{'assembly_span'};
	my $atgc = $stats->{'assembly_atgc'};
	my $gc = $stats->{'assembly_gc'} * 100;
	my $n = $span - $atgc;
	my $cegma_complete = $stats->{'cegma_complete'};
	my $cegma_partial = $stats->{'cegma_partial'};
	my @lengths = @{$stats->{'lengths'}};
	my @ctg_lengths = @{$stats->{'contig_lengths'}};
	my $scaffolds = 'scaffolds: ['.join(',',@lengths).']';
	my $contigs = @ctg_lengths ? 'contigs: ['.join(',',@ctg_lengths).'],' : '';
	my $gcs = $stats->{'GCs'} ? 'GCs: ['.join(',',@{$stats->{'GCs'}}).'],' : '';
	my $ns = $stats->{'Ns'} ? 'Ns: ['.join(',',@{$stats->{'Ns'}}).'],' : '';
	open ASM,">summary/$production_name"."_assembly.html";
print ASM <<EOF;
<script type="text/javascript" src="/components/00_jquery.min.js"></script>
<script type="text/javascript" src="/components/00_d3.js"></script>
<script type="text/javascript" src="/components/99_assembly_stats.js"></script>
<script>
  var stats = { assembly: $span,
                ATGC:     $atgc,
			    GC:       $gc,
			    N:		  $n,
			    $gcs
			    $ns
			    cegma_complete:	$cegma_complete,
			    cegma_partial:	$cegma_partial,
			    $contigs
			    $scaffolds
			};
  var asm = new Assembly (stats);
  \$(document).ready(function(){
      var svg = d3.select('#assembly_stats')
        		.append('svg');

      asm.drawPlot(svg);
    });
</script>
<div id="assembly_stats"></div>

EOF

	close ASM;
	return 1;
}


sub web_ini_file {
	my $params = shift;
	my $stats = shift;
	my $production_name = $params->{'META'}{'SPECIES.PRODUCTION_NAME'};
	my $species = $params->{'META'}{'SPECIES.SCIENTIFIC_NAME'};
	my $filename = $params->{'DATABASE_CORE'}{'NAME'};
	my $version = $params->{'META'}{'GENEBUILD.VERSION'} || 1;
	my $size = int($stats->{'assembly_span'}/10000000)/100 ;
	my $location_param = $params->{'META'}{'SAMPLE.LOCATION_PARAM'} || '';
	my $location_text = $params->{'META'}{'SAMPLE.LOCATION_TEXT'} || '';
	my $gene_param = $params->{'META'}{'SAMPLE.GENE_PARAM'} || '';
	my $gene_text = $params->{'META'}{'SAMPLE.GENE_TEXT'} || '';
	my $transcript_param = $params->{'META'}{'SAMPLE.TRANSCRIPT_PARAM'} || '';
	my $transcript_text = $params->{'META'}{'SAMPLE.TRANSCRIPT_TEXT'} || '';
	my $search_text = $params->{'META'}{'SAMPLE.SEARCH_TEXT'} || 'Peptidoglycan';
	open INI,">summary/$production_name.ini";
print INI <<EOF;
[general]
SPECIES_RELEASE_VERSION = $version
ENSEMBL_GENOME_SIZE     = $size

[databases]
DATABASE_CORE = $filename
DATABASE_USERDATA  = userdata

[DATABASE_USERDATA]

[ENSEMBL_STYLE]

[ENSEMBL_EXTERNAL_URLS]

[ENSEMBL_SPECIES_SITE]

[SPECIES_DISPLAY_NAME]

[ENSEMBL_EXTERNAL_DATABASES]

[ENSEMBL_EXTERNAL_INDEXERS]

[S4_PROTEIN]

[S4_PROTEIN_STRUCTURE]

[S4_LITERATURE]

[SAMPLE_DATA]
LOCATION_PARAM    = $location_param
LOCATION_TEXT     = $location_text

GENE_PARAM        = $gene_param
GENE_TEXT         = $gene_text

TRANSCRIPT_PARAM  = $transcript_param
TRANSCRIPT_TEXT   = $transcript_text

SEARCH_TEXT       = $search_text
EOF

close INI;
	return 1;
}

sub readable {
	my ($number,$unit) = @_;
	my $exponent = 0;
	while ($number / 10 > 1){
		$number /= 10;
		$exponent++;
	}
	my %suffix = ( 	0	=> ''.$unit,
					3	=> 'k'.$unit,
					6	=> 'M'.$unit,
					9	=> 'G'.$unit,
					12	=> 'T'.$unit,
					15	=> 'P'.$unit);
	$number = sprintf("%.2f", $number);
	while ($exponent % 3 > 0){
		$number *= 10;
		$exponent--;
	}
	return "$number $suffix{$exponent}";
}

sub gff_feature_summary {
	my ($params,$infile) = @_;
	my %features_by_type;
	my $file = $infile;
	if ($params->{'GFF'}{'SORT'} && $file !~ m/\.sorted/){
		$file .= ".sorted";
		if (!-e $file){
			system "cat ".$infile." | sort > ".$file.".tmp";
			system "mv ".$file.".tmp ".$file;
		}
	}
	unshift @ARGV,$file;
	my $gff = GFFTree->new({});
	$gff->name('GFF');
	$gff->expect_columns(9,'skip');
	$gff->lacks_id('all','make');
	$gff->lacks_id('cds','Name');
	my @tsc_types;
	foreach my $key (keys %{$params->{'GFF'}}){
		my $value = $params->{'GFF'}{$key};
		if (ref $value || ref $value eq 'ARRAY') {
			my @value = @$value;
			my $type = shift @value;
			if ($type eq 'MULTILINE'){
				$gff->multiline(@value);
			}
			elsif ($type eq 'EXPECTATION'){
				if ($value[1] =~ m/^hasParent$/i && $value[2] =~ m/^gene$/i){
					push @tsc_types, $value[0];
				}
			}
		}
	}
	while ($gff->parse_chunk(@{$params->{'GFF'}{'CHUNK'}})){
		my @features = $gff->daughters();
		foreach my $feature (@features){
			_summarise_feature($feature,\%features_by_type);
		}
	}
	# print summaries to files
	mkdir "summary";
	open ATTR,">summary/GFF.$file.attribute_counts.txt";
	foreach my $type (sort keys %features_by_type){
		foreach my $attribute (sort keys %{$features_by_type{$type}{'attributes'}}){
			print ATTR "$type\t$attribute\t$features_by_type{$type}{'attributes'}{$attribute}\n";
		}
	}
	close ATTR;

	my %return_stats;
	foreach my $type (sort keys %features_by_type){
		if ($type =~ m/^gene$/i){
			$return_stats{'gene_count'} = $features_by_type{$type}{'attributes'}{'ID'};
		}
		elsif ($type =~ m/^cds$/i){
			$return_stats{'translation_count'} = $features_by_type{$type}{'attributes'}{'ID'};
		}
		else {
			foreach my $tsc (@tsc_types){
				if ($type =~ m/^$tsc$/i){
					$return_stats{'transcript_count'} += $features_by_type{$type}{'attributes'}{'ID'};
				}
			}
		}
	}
	return \%return_stats;
}


sub prepared_gff_feature_summary {
	my ($params,$infile,$features) = @_;
	my %features_by_type = %$features if $features;
	my $production_name = $params->{'META'}{'SPECIES.PRODUCTION_NAME'};
	my $file = $infile;
	unshift @ARGV,$file;
	my $gff = GFFTree->new({});
	$gff->name('GFF');
	$gff->expect_columns(9,'skip');
	$gff->lacks_id('all','make');
	$gff->lacks_id('cds','Name');
	my @tsc_types;
	foreach my $key (keys %{$params->{'GFF'}}){
		my $value = $params->{'GFF'}{$key};
		if (ref $value || ref $value eq 'ARRAY') {
			my @value = @$value;
			my $type = shift @value;
			if ($type eq 'MULTILINE'){
				$gff->multiline(@value);
			}
			elsif ($type eq 'EXPECTATION'){
				if ($value[1] =~ m/^hasParent$/i && $value[2] =~ m/^gene$/i){
					push @tsc_types, $value[0];
				}
			}
		}
	}
	while ($gff->parse_chunk(@{$params->{'GFF'}{'CHUNK'}})){
		my @features = $gff->daughters();
		foreach my $feature (@features){
			_summarise_feature($feature,\%features_by_type);
		}
	}
	# print summaries to files
	mkdir "summary";
	open ATTR,">summary/GFF.$file.attribute_counts.txt";
	foreach my $type (sort keys %features_by_type){
		foreach my $attribute (sort keys %{$features_by_type{$type}{'attributes'}}){
			print ATTR "$type\t$attribute\t$features_by_type{$type}{'attributes'}{$attribute}\n";
		}
	}
	close ATTR;
	open RDATA,">summary/${production_name}_feature";
	foreach my $type (sort keys %features_by_type){
		for (my $i = 0; $i < @{$features_by_type{$type}{'lengths'}}; $i++){
			print RDATA "$type\t$features_by_type{$type}{'regions'}[$i]\t$features_by_type{$type}{'sources'}[$i]\t$features_by_type{$type}{'lengths'}[$i]\n";
		}
	}
	close RDATA;
	make_feature_plot("summary/${production_name}_feature");

	my %return_stats;
	foreach my $type (sort keys %features_by_type){
		if ($type =~ m/^gene$/i){
			$return_stats{'gene_count'} = $features_by_type{$type}{'attributes'}{'ID'};
		}
		elsif ($type =~ m/^cds$/i){
			$return_stats{'translation_count'} = $features_by_type{$type}{'attributes'}{'ID'};
		}
		else {
			foreach my $tsc (@tsc_types){
				if ($type =~ m/^$tsc$/i){
					$return_stats{'transcript_count'} += $features_by_type{$type}{'attributes'}{'ID'};
				}
			}
		}
	}
	return (\%return_stats,\%features_by_type);
}


sub cegma_file_summary {
	my ($infile) = @_;
	my %return_stats;
	open CEGMA,"$infile";
	while (<CEGMA>){
		chomp;
		if (m/^\s*Complete\s+\S+\s+(\S+)\s/){
			$return_stats{'cegma_complete'} = $1;
		}
		elsif (m/^\s*Partial\s+\S+\s+(\S+)\s/){
			$return_stats{'cegma_partial'} = $1;
		}
	}
	close CEGMA;
	return (\%return_stats);
}



sub make_feature_plot {

my ($rdatafile) = @_;
open RSCRIPT, ">$rdatafile.R" or die $!;

print RSCRIPT '
args <- commandArgs(trailingOnly = TRUE)
data=read.delim(args[1],header=F,col.names=c("type","region","source","len"))
distrib=aggregate(data$len, by=list(data$region,data$type), FUN=length);
x = c(1,2,3,4,5)
j = 6
names(x) = c("gene","mRNA","exon","CDS")
cols = c(rgb(31/255,120/255,180/255),rgb(166/255,206/255,227/255),rgb(51/255,160/255,44/255),
rgb(251/255,154/255,153/255),rgb(178/255,223/255,138/255),rgb(227/255,26/255,28/255),
rgb(253/255,191/255,111/255),rgb(255/255,127/255,0/255),rgb(202/255,178/255,214/255))
alltypes = unique(unlist(data[,"type"], use.names = FALSE))
allsources = unique(unlist(data[,"source"], use.names = FALSE))
types = c(NA)
colormap = c(NA)
for (i in 1:length(alltypes)){
	if (sum(data[,"type"] == alltypes[i]) > 10){
		types = c(as.character(alltypes[i]),types)
		if (!is.na(x[types[1]])){
			colormap = c(x[types[1]],colormap)
		}
		else {
			colormap = c(j,colormap)
			j = j+1
		}
	}
}
types <- c("gene","mRNA","exon","CDS")

height=150+150*length(types)
png(paste(args[1],"dist","png",sep="."),900,height,res=200)
distrib[distrib[,3]<1,3] <- 0.9
par(las=1,mai=c(0,0,0.1,0),oma=c(4,4,1,1),cex.axis=0.8,cex.lab=1,mfrow=c(4,1))

for (i in 1:length(types)){
	h = hist(log(distrib[distrib[,2]==types[i],3],base=10),breaks=seq(from=0,to=5,by=0.25),plot=F)
	h$density = h$counts/sum(h$counts)*100
	h$density[h$density>0 & h$density<1] = 1
	plot(h,col=cols[i],border="white",freq=F,ylab="",xlab="",main="",axes=F)
	axis(2)
	legend("bottomright",c(types[i],paste("n = ",sum(h$counts))),ncol=1,box.col=NA,bg=NA)
	if (i == 1){
		#mtext(text="GFF feature distributions by region", side=3, line=0, las=0, cex=1)
	}
	if (i == ceiling(length(types)/2)){
		mtext(text="Percentage", side=2, line=2.5, las=0, cex=0.8)
	}
	if (i == length(types)){
		axis(1,labels=c("1e+00","1e+01","1e+02","1e+03","1e+04","1e+05"),at=seq(from=0,to=5,by=1))
		mtext(text="Count per region", side=1, line=2.5, las=0, cex=0.8)
	}
	else {
		axis(1,labels=F)
	}
}
dev.off()
png(paste(args[1],"length","png",sep="."),900,height,res=200)
data[data[,"len"]<1,"len"] <- 0.9
par(las=1,mai=c(0,0,0.1,0),oma=c(4,4,1,1),cex.axis=0.8,cex.lab=1,mfrow=c(4,1))

for (i in 1:length(types)){
	h = hist(log(data[data$type==types[i],"len"],base=10),breaks=seq(from=-0.5,to=6,by=0.5),plot=F)
	h$density = h$counts/sum(h$counts)*100
	h$density[h$density>0 & h$density<1] = 1
	plot(h,col=cols[i],border="white",freq=F,ylab="",xlab="",main="",axes=F)
	axis(2)
	legend("bottomleft",c(types[i],paste("n = ",sum(h$counts))),ncol=1,box.col=NA,bg=NA)
	if (i == 1){
		#mtext(text="GFF feature length distributions", side=3, line=0, las=0, cex=1)
	}
	if (i == ceiling(length(types)/2)){
		mtext(text="Percentage", side=2, line=2.5, las=0, cex=0.8)
	}
	if (i == length(types)){
		axis(1,labels=c("1e+00","1e+01","1e+02","1e+03","1e+04","1e+05","1e+06"),at=seq(from=0,to=6,by=1))
		mtext(text="Length (bp)", side=1, line=2.5, las=0, cex=0.8)
	}
	else {
		axis(1,labels=F)
	}
}
dev.off()
';
close RSCRIPT;
system "Rscript $rdatafile.R $rdatafile >/dev/null";
#unlink  "$rdatafile.R";
}

sub _summarise_feature {
	my ($feature,$features_by_type) = @_;
	my $region = $feature->_seq_name() || 'NA';
	my $source = $feature->{attributes}->{_source} || 'NA';
	my $type = $feature->{attributes}->{_type} || 'NA';
	my $length = $feature->_length  || 0;
	push @{$features_by_type->{$type}{'lengths'}},$length;
	push @{$features_by_type->{$type}{'sources'}},$source;
	push @{$features_by_type->{$type}{'regions'}},$region;
	foreach my $attribute (keys %{$feature->{attributes}}){
		next if $attribute =~ m/^_/;
		next if $attribute =~ m/_array$/;
		$features_by_type->{$type}{'attributes'}{$attribute}++;
	}

	my @features = $feature->daughters();
	foreach my $subfeature (@features){
		_summarise_feature($subfeature,$features_by_type);
	}
	return 1;
}

sub split_gff {
	my ($params,$infiles) = @_;
	my $pattern = $params->{'GFF'}{'SPLIT'}[0];
	system "csplit ".$infiles->{$params->{'GFF'}{'SPLIT'}[1]}{'name'}." /$pattern/";
	if ($params->{'GFF'}{'SPLIT'}[0] && -e "xx00"){
		system "mv xx00 ".$infiles->{$params->{'GFF'}{'SPLIT'}[1]}{'name'};
	}
	else {
		die "ERROR: [GFF] SPLIT failed to generate a file of type ".$infiles->{$params->{'GFF'}{'SPLIT'}[1]}{'type'}."\n";
	}
	if ($params->{'GFF'}{'SPLIT'}[1] && -e "xx01"){
		system "mv xx01 ".$infiles->{$params->{'GFF'}{'SPLIT'}[2]}{'name'};
	}
	else {
		die "ERROR: [GFF] SPLIT failed to generate a file of type ".$infiles->{$params->{'GFF'}{'SPLIT'}[2]}{'type'}."\n";
	}

	return 1;
}

sub fetch_file {
	my ($type,$location,$new_name);
	$location = shift;
	if (ref $location){
		$type = $location->[0];
		if (defined $location->[2]){
			$new_name = $location->[2];
		}
		$location = $location->[1];
	}
	# work out file name from location
	$location =~ m/.+\/([^\/]+)$/;
	my $filename = $1 ? $1 : $location;
	my $command;
  my $compression;
	if ($filename =~ s/\.(gz|gzip|tar\.gz|tgz|zip)$//){
		$compression = ".".$1;
	}
	if (($new_name && !-e $new_name) || (!$new_name && !-e $filename.$compression)){
		if ($location =~ m/^(?:ftp|http|https):/){
			$command = 'wget';
			system "wget $location -O $filename"."$compression";
		}
		elsif ($location =~ m/:[\/~]/){
			$command = 'scp';
			system "scp $location $filename"."$compression";
		}
		else {
			$command = 'cp';
			system "cp $location $filename"."$compression";
		}
	}
	if ($compression && !-e $filename){
		if (!-e $filename.$compression){
			# something went wrong when copying the file
			die "ERROR: could not $command $location to $filename"."$compression\n";
		}
		else {
			if ($compression =~ m/^\.t/){
				system "tar xf $filename"."$compression";
			}
			elsif ($compression =~ m/^\.g/){
				system "gunzip $filename"."$compression";
			}
			elsif ($compression =~ m/^\.z/){
				system "unzip $filename"."$compression";
			}
			if (!-e $filename){
				# this compression type is not currently supported
				die "ERROR: could not extract $filename"."$compression to $filename\n";
			}
		}
	}
	if (!-e $filename){
		# something went wrong when copying the file
		die "ERROR: could not $command $location to $filename\n";
	}
	if (!$type){
		# TODO: infer type properly
		if ($filename =~ m/\.([^\.])$/){
			$type = $1;
			warn "WARNING: no type specified for $filename from $location, inferring '$type' based on file extension\n";

		}
		else {
			die "ERROR: no type specified for $filename from $location, unable to infer based on file extension\n";
		}
	}
	if (defined $new_name){
		system "mv $filename $new_name";
		if (!-e $new_name){
			die "ERROR: could not move $filename from $location to $new_name\n";
		}
		$filename = $new_name;
	}
	if ($type =~ m/^(?:fa|faa|fas|fasta|fna|fsa)/i  ){
		$type = 'fas';
	}
	elsif ($type =~ m/^gff/i){
		$type = 'gff';
	}
	elsif ($type =~ m/^agp/i){
		$type = 'agp';
	}
	elsif ($type =~ m/^csv/i){
		$type = 'csv';
	}
	elsif ($type =~ m/^tsv/i){
		$type = 'tsv';
	}
	return ($filename,$type);
}

sub template_db_connect {
	my $params = shift;
	my $dsn = "DBI:mysql:host=$params->{'DATABASE_TEMPLATE'}{'HOST'};port=$params->{'DATABASE_TEMPLATE'}{'PORT'}";
	my $dbh = DBI->connect($dsn,"$params->{'DATABASE_TEMPLATE'}{'RO_USER'}") || die "ERROR: unable to connect to [DATABASE_TEMPLATE] HOST ".$params->{'DATABASE_TEMPLATE'}{'HOST'}." using provided settings";
	$dbh->do('USE '.$params->{'DATABASE_TEMPLATE'}{'NAME'}) || die "ERROR: unable to connect to [DATABASE_TEMPLATE] NAME ".$params->{'DATABASE_TEMPLATE'}{'NAME'}." using provided settings";
	return $dbh;
}

sub taxonomy_db_connect {
	my $params = shift;
	my $dsn = "DBI:mysql:host=$params->{'DATABASE_TAXONOMY'}{'HOST'};port=$params->{'DATABASE_TAXONOMY'}{'PORT'}";
	my $dbh = DBI->connect($dsn,"$params->{'DATABASE_TAXONOMY'}{'RO_USER'}") || die "ERROR: unable to connect to [DATABASE_TAXONOMY] HOST ".$params->{'DATABASE_TAXONOMY'}{'HOST'}." using provided settings";
	$dbh->do('USE '.$params->{'DATABASE_TAXONOMY'}{'NAME'}) || die "ERROR: unable to connect to [DATABASE_TAXONOMY] NAME ".$params->{'DATABASE_TAXONOMY'}{'NAME'}." using provided settings";
	return $dbh;
}

sub core_db_host_connect {
	my $params = shift;
	my $dsn = "DBI:mysql:host=$params->{'DATABASE_CORE'}{'HOST'};port=$params->{'DATABASE_CORE'}{'PORT'}";
	my $dbh = DBI->connect($dsn,"$params->{'DATABASE_CORE'}{'RW_USER'}","$params->{'DATABASE_CORE'}{'RW_PASS'}") || die "ERROR: unable to connect to [DATABASE_CORE] HOST ".$params->{'DATABASE_CORE'}{'HOST'}." using provided settings";
	return $dbh;
}

sub core_db_connect {
	my $params = shift;
	my $dbh = core_db_host_connect($params);
	$dbh->do('CREATE DATABASE IF NOT EXISTS '.$params->{'DATABASE_CORE'}{'NAME'}) || die "ERROR: unable to create [DATABASE_CORE] NAME ".$params->{'DATABASE_CORE'}{'NAME'}." using provided settings";
	$dbh->do('USE '.$params->{'DATABASE_CORE'}{'NAME'}) || die "ERROR: unable to connect to [DATABASE_CORE] NAME ".$params->{'DATABASE_CORE'}{'NAME'}." using provided settings";
	return $dbh;
}


sub clone_ensembl_code {
	my $params = shift;
	if (!-d $params->{'ENSEMBL'}{'LOCAL'}){
		system 'mkdir -p '.$params->{'ENSEMBL'}{'LOCAL'};
		if (!-d $params->{'ENSEMBL'}{'LOCAL'}){
			die 'ERROR: unable to create directory [ENSEMBL] LOCAL '.$params->{'ENSEMBL'}{'LOCAL'}."\n";
		}
	}
	if (!-d $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl'){
		system 'git clone -b '.$params->{'ENSEMBL'}{'BRANCH'}.' '.$params->{'ENSEMBL'}{'ENSEMBL'}.' '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl';
		if (!-d $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl'){
			die 'ERROR: unable to checkout [ENSEMBL] code BRANCH '.$params->{'ENSEMBL'}{'BRANCH'}.' from ENSEMBL repository '.$params->{'ENSEMBL'}{'ENSEMBL'}."\n";
		}
	}
	if (!-d $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-production'){
		system 'git clone -b '.$params->{'ENSEMBL'}{'BRANCH'}.' '.$params->{'ENSEMBL'}{'PRODUCTION'}.' '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-production';
		if (!-d $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-production'){
			die 'ERROR: unable to checkout [ENSEMBL] code BRANCH '.$params->{'ENSEMBL'}{'BRANCH'}.' from PRODUCTION repository '.$params->{'ENSEMBL'}{'PRODUCTION'}."\n";
		}
	}
	if (!-d $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline'){
		system 'git clone -b master '.$params->{'ENSEMBL'}{'PIPELINE'}.' '.$params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline';
		if (!-d $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline'){
			die 'ERROR: unable to checkout [ENSEMBL] code from PIPELINE repository '.$params->{'ENSEMBL'}{'PIPELINE'}."\n";
		}
	}
	if (!-e $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl/sql/table.sql'){
		die "ERROR: unable to find an sql file in the [ENSEMBL] ENSEMBL LOCAL directory ".$params->{'ENSEMBL'}{'LOCAL'}."/ensembl/sql/\n";
	}
	if (!-e $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-pipeline/sql/table.sql'){
		die "ERROR: unable to find an sql file in the [ENSEMBL] PIPELINE LOCAL directory ".$params->{'ENSEMBL'}{'LOCAL'}."/ensembl-pipeline/sql/\n";
	}
	return 1;
}

sub load_ini {
	my ($params,$ini_file,$sections) = @_;
	open INI,$ini_file || die "ERROR: unable to open file $ini_file\n",usage(),"\n";
	my %sections = %$sections;
	my $section;
	while (<INI>){
		chomp;
		# treat anything after a semicolon as a comment
		s/;.+//;
		# unescape semicolons
		s/\%3B/;/g;
		if (/^\s*\[(.+)\]\s*$/){
			$section = $1;
		}
		elsif (/^\s*(.+?)\s*=\s*(.+?)\s*$/) {
			my $key = $1;
			my $value = $2;
			if ($value =~ m/^\[\s*(.+)\s*\]/){
				my @value = split /\s+/,$1;
				# unescape spaces
				for (my $i = 0; $i < @value; $i++){
					$value[$i] =~ s/\%20/ /g;
				}
				$value = \@value;
			}
			$params->{$section}{$key} = $value;
		}
	}
	my $warnings = 0;
	foreach $section (keys %sections){
		if (scalar keys %{$sections{$section}} > 0){
			if (!$params->{$section}){
				# warn missing section
				warn "WARNING: missing section [$section]\n";
				$warnings++;
			}
			foreach my $subsection (keys %{$sections{$section}}){
				if ($sections{$section}{$subsection} == 1){
					if (!$params->{$section}{$subsection}){
						# warn missing subsection
						warn "WARNING: missing subsection [$section] $subsection\n";
						$warnings++;
					}
				}
				else {
					if (!$params->{$section}{$subsection}){
						# look for alternate section
						my $alt;
						for (my $i = 0; $i < @{$sections{$section}{$subsection}}; $i++){
							if ($params->{$section}{$sections{$section}{$subsection}->[$i]}){
								$alt = 1;
								last;
							}
						}
						next if $alt;
						my $str = join (' or ',@{$sections{$section}{$subsection}});
						warn "WARNING: section [$section] must have at least one subsection of type $subsection or $str\n";
						$warnings++;
					}
				}
			}
		}
	}
	# now check that all files in GENE_DESCRIPTIONS, GENE_NAMES, TRANSCRIPT_DESCRIPTIONS and TRANSCRIPT_NAMES also have GENE_STABLE_IDS/TRANSCRIPT_STABLE_IDS and FILES
	my @features = qw ( GENE TRANSCRIPT );
	my @properties = qw ( DESCRIPTIONS NAMES );
	foreach my $feature (@features){
		foreach my $property (@properties){
			foreach my $subsection (keys %{$params->{$feature.'_'.$property}}){
				if (!$params->{$feature.'_STABLE_IDS'}{$subsection}){
					# warn no stable_id for property
					warn "WARNING: no matching [".$feature."_STABLE_IDS] $subsection for [".$feature."_".$property."] $subsection\n";
					$warnings++;
				}
				if (!$params->{'FILES'}{$subsection}){
					# warn no file for property
					warn "WARNING: no matching [FILES] $subsection for [".$feature."_".$property."] $subsection\n";
					$warnings++;
				}
			}
		}
	}
	if ($warnings > 0){
		die "ERROR: unable to parse ini file $ini_file without warnings\n";
	}

	return 1;
}

sub table_exists {
	#http://www.perlmonks.org/?node_id=284436#tablecheck
    my $db = shift;
    my $table = shift;
    my @tables = $db->tables('','','','TABLE');
    if (@tables) {
        for (@tables) {
            next unless $_;
            return 1 if $_ =~ m/$table/;
        }
    }
    else {
        eval {
            local $db->{PrintError} = 0;
            local $db->{RaiseError} = 1;
            $db->do(qq{SELECT * FROM $table WHERE 1 = 0 });
        };
        return 1 unless $@;
    }
    return 0;
}



1;
