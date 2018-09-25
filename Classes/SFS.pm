package SFS;

use strict;
use warnings;
use Classes::Misc;


sub new{
	my $class=shift;
	my $self={};
	bless $self, $class;
	return $self;
	}


sub printSNPDifferentials{
	my ($self, $outfile, $regions, $ref_annotation, $refseq, $ancref, $chrommap, $ref_genotypes, $ref_popsizes, $ref_groups, $ref_mingroupind, $mincov, $haplotize, $ref_regionstomask, $onlyannotated)=@_;

	die "Script only supports a single pair of groups!\n" unless (@$ref_groups==2);
	open my $output, ">>", $outfile or die "Could not open output file $outfile!$!\n";
	my @type_of_change=('W2W', 'W2S', 'S2W', 'S2S');

	ID: for my $id ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($id) ){
			my $start=$regions->{_Loci}{$id}[$locus]{'start'};
			my $end=$regions->{_Loci}{$id}[$locus]{'end'};
			my ($startchrom, $endchrom, $trstart, $trend);
			if (defined $chrommap){ 
				($startchrom, $trstart)=$chrommap->findPosonChrom($id, $start);
				($endchrom, $trend)=$chrommap->findPosonChrom($id, $end);
				if (defined $startchrom && defined $endchrom && $startchrom ne $endchrom){ die "Start and end position of window $id:$start-$end are not on the same chromosome!\n" }
				}
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($id, $start, $end, 1);
			my $ref_ancstring=$ancref->randomAccess($id, $start, $end, 1);
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($id) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($id) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($hardmasked, $notbiallelic, $undefanc, $undefcat, $filtered)=(0, 0, 0, 0, 0);
			my $nsnps_tot=0;
			my $prevpos=$start;

			POSITION: for my $pos ( @$ref_range ){
				my $gap=$pos-$prevpos;
				substr($$ref_seqstring, 0, $gap, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				substr($$ref_ancstring, 0, $gap, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				my $refbase=uc( substr($$ref_seqstring, 0, 1) );
				my $ancbase=uc( substr($$ref_ancstring, 0, 1) );
				if ($gap>1){ $hardmasked+=$gap-1 }
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				if ($ancbase eq 'N' ){ ++$undefanc; next POSITION }	# discard if base in ancestral reference sequence is undefined.

				next POSITION unless ( defined $ref_genotypes->getValue($id, $pos) );	# only consider polymorphic sites.
				next POSITION if ($onlyannotated && !exists $ref_annotation->{_Annotations}{$id}{$pos});	# only consider exonic sites.
				++$nsnps_tot;

				my $ancallele=$ref_genotypes->getAlleleforBase($id, $pos, $ancbase);
				my ($ancestral, @derived)=$ref_genotypes->getPolarizedAlleles($id, $pos, $ancallele);
				unless (defined $ancestral){ warn "$id:$pos - no ancestral allele for ancestral base $ancbase!\n"; ++$undefanc; next POSITION }
				$derived[0]=!$ancestral unless (@derived);
				if (@derived>1){ ++$notbiallelic; next POSITION }	# maximum one derived allele permitted.
				my $derbase=$ref_genotypes->getBase($id, $pos, $derived[0]);
				my $cat=Misc::getWeakStrong($ancbase, $derbase);
				unless (defined $cat){ ++$undefcat; next POSITION }

				my @nchrom=(0) x @$ref_groups;
				my @nder=(0) x @$ref_groups;

				GROUP: for my $group ( 0..$#{$ref_groups} ){
					my @validalleles;
					for my $subpop ( @{$ref_groups->[$group]} ){
						IND: for my $ind (0..$ref_popsizes->[$subpop]-1){
							my $coverage=$ref_genotypes->getCoverage($id, $pos, $subpop, $ind);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							if ($haplotize){
								my $index=2*$ind+int(rand(2));
								my $allele=$ref_genotypes->getValue($id, $pos, $subpop, $index);
								unless (defined $allele && $allele ne '.'){ next IND }
								push @validalleles, $allele;
								}
							else {
								my $index=2*$ind;
								my $allele1=$ref_genotypes->getValue($id, $pos, $subpop, $index);
								unless (defined $allele1 && $allele1 ne '.'){ next IND }
								++$index;
								my $allele2=$ref_genotypes->getValue($id, $pos, $subpop, $index);
								unless (defined $allele2 && $allele2 ne '.'){ next IND }
								push @validalleles, ($allele1, $allele2);
								}
							}
						}

					my $minnumber=$haplotize ? $ref_mingroupind->[$group] : 2*$ref_mingroupind->[$group];
					unless (@validalleles>=$minnumber){ ++$filtered; next POSITION }
					for my $allele (@validalleles){ ++$nder[$group] if ($allele ne $ancestral) }
					$nchrom[$group]=@validalleles;
					}

				my ($trchrom, $trpos);
				if (defined $chrommap){
					($trchrom, $trpos)=$chrommap->findPosonChrom($id, $pos);
					}
				++$trpos if (defined $trpos);	# change to 1-based indexing.
				my $higher_group=($nder[0]/$nchrom[0]>$nder[1]/$nchrom[1]) ? 0 : 1;
				$higher_group="-" if ($nder[0]/$nchrom[0]==$nder[1]/$nchrom[1]);
				print $output "$id\t", $pos+1, "\t", $trchrom // 'NA', "\t", $trpos // 'NA',
					"\t$ancbase\t$derbase\t$type_of_change[$cat]\t", $ref_annotation->{_Annotations}{$id}{$pos} // 'noncoding',
					"\t", sprintf("%.4f", abs($nder[0]/$nchrom[0]-$nder[1]/$nchrom[1]) ), "\t$higher_group";
				for my $group ( 0..$#{$ref_groups} ){
					print $output "\t", sprintf("%.4f", $nder[$group]/$nchrom[$group]), "\t$nder[$group]\t", $nchrom[$group]/2;
					}
				print $output "\n";	# scaffold    position    chromosome    position    type_of_change     allele_frequency_difference   group_index_with_higher_derived_frequency {derived_allele_frequency    derived_allele_count   number_of_individuals_covered}*2
				} continue { $prevpos=$pos }
			print "$id:$start-$end; total number of SNPs: $nsnps_tot, hardmasked: $hardmasked, not biallelic: $notbiallelic, undefined ancestral state: $undefanc, undefined category: $undefcat, filtered: $filtered\n";
			}
		}
	close $output;
	return;
	}


sub validRegions{
	my ($self, $outfile, $regions, $refseq, $ancref, $ref_missing, $ref_popsizes, $ref_minpopind, $ref_groups, $ref_mingroupind, $mincov, $verbose)=@_;

	my $ref_popsums=$ref_missing->getPopsums();

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">>", $outfile or die "Could not open output file!\n" }
	
	my $validtot=0;

	ID: for my $id ( $regions->allKeys() ){
		LOCUS: for my $locus ( sort { $a->{'start'}<=>$b->{'start'} } $regions->allValues($id) ){
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($id, $start, $end, 1);
			my $ref_ancstring=(defined $ancref) ? $ancref->randomAccess($id, $start, $end, 1) : \$$ref_seqstring;
			my ($totsites, $validsites, $hardmasked, $undefanc, $insuffcov)=(0, 0, 0, 0, 0);
			my ($startpos, $endpos, $valid);

			POSITION: for my $pos ($start..$end){
				++$totsites;
				$valid=0;
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );
				my $ancbase=uc( substr($$ref_ancstring, 0, 1, "") );
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				if ($ancbase eq 'N' ){ ++$undefanc; next POSITION }	# discard if base in ancestral reference sequence is undefined.
				my $offset=$pos-$start;

				GROUP: for my $group ( 0..$#{$ref_groups} ){
				    my $validgroupind=0;
					POP: for my $subpop ( @{$ref_groups->[$group]} ){
    					my $validpopind=0;
    					IND: for my $ind (0..$ref_popsizes->[$subpop]-1){
    						my $allind=$ind+$ref_popsums->[$subpop];
							my $coverage=$ref_missing->getValue($id, $start, $allind, $offset);
    						unless (defined $coverage && $coverage>=$mincov){ next IND }							
    						++$validpopind;
    						}
    					unless ($validpopind>=$ref_minpopind->[$subpop]){ ++$insuffcov; next POSITION }
    					$validgroupind+=$validpopind;
    				    }
   					unless ($validgroupind>=$ref_mingroupind->[$group]){ ++$insuffcov; next POSITION }
    				}
    			$valid=1;
    		    ++$validsites;
    		    $startpos=$pos unless (defined $startpos);
    		    $endpos=$pos;
                } continue { unless ($valid){ print $output "$id\t$startpos\t", $endpos+1, "\n" if (defined $startpos); $startpos=undef } }
			print $output "$id\t$startpos\t", $endpos+1, "\n" if (defined $startpos);
			$validtot+=$validsites;
			print STDERR "$id:$start-$end, total sites considered: $totsites, valid sites: $validsites, hardmasked: $hardmasked, undefined ancestral state: $undefanc, insufficent coverage: $insuffcov\n" if $verbose;
            }
        }
	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return($validtot);
	}


1;

