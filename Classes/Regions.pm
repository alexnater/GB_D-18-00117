package Regions;

use strict;
use warnings;
use Classes::Misc;
use List::Util qw( min );


sub new{
	my $class=shift;
	my $loci=shift;
	my $self={};
	$self->{_Loci}=defined $loci ? $loci : {};
	bless $self, $class;
	return $self;
	}


sub getValue{
	my ($self, $chrom, $locus, $key)=@_;

	if (@_==1){return($self->{_Loci})}
	elsif (@_==2){return($self->{_Loci}{$chrom})}
	elsif (@_==3){return($self->{_Loci}{$chrom}[$locus])}
	elsif (@_==4){return($self->{_Loci}{$chrom}[$locus]{$key})}
	else {die "Invalid number of arguments!\n"}
	}


sub setValue{
	my $self=shift;
	unless (@_==4){die "Wrong number of arguments!\n"}
	my ($chrom, $locus, $key, $value)=@_;

	$self->{_Loci}{$chrom}[$locus]{$key}=$value;

	return($self->{_Loci}{$chrom}[$locus]{$key});
	}


sub getExcluded{
	my $self=shift;
	my ($chrom, $locus)=@_;

	if (@_==2){ return( $self->{_Loci}{$chrom}[$locus]{'excluded'} ) }
	else { die "Invalid number of arguments!\n" }
	}


sub setExcluded{
	my $self=shift;
	unless (@_==3){ die "Wrong number of arguments!\n" }
	my ($chrom, $locus, $value)=@_;
	unless ($value==0 || $value==1){ die "Invalid value, 0 or 1 required!\n" }

	$self->{_Loci}{$chrom}[$locus]{'excluded'}=$value;

	return( $self->{_Loci}{$chrom}[$locus]{'excluded'} );
	}


sub getSize{
	my ($self, $chrom, $locus)=@_;

	if (@_==1){return(scalar(keys %{$self->{_Loci}}))}
	elsif (@_==2){return(scalar @{$self->{_Loci}{$chrom}})}
	elsif (@_==3){return(scalar(keys %{$self->{_Loci}{$chrom}[$locus]}))}
	else {die "Invalid number of arguments!\n"}
	}


sub allKeys{
	my ($self, $chrom, $locus)=@_;

	if (@_==1){return(keys %{$self->{_Loci}})}
	elsif (@_==2){return((0..$#{$self->{_Loci}{$chrom}}))}
	elsif (@_==3){return(keys %{$self->{_Loci}{$chrom}[$locus]})}
	else {die "Invalid number of arguments!\n"}
	}


sub allValues{
	my ($self, $chrom, $locus)=@_;

	if (@_==1){return(values %{$self->{_Loci}})}
	elsif (@_==2){return(@{$self->{_Loci}{$chrom}})}
	elsif (@_==3){return(values %{$self->{_Loci}{$chrom}[$locus]})}
	else {die "Invalid number of arguments!\n"}
	}


sub findLocus{
	my ($self, $id, $pos)=@_;
	unless (@_==3){die "Wrong number of arguments!\n"}

	my $index;
	LOCUS: for my $locus (sort {$self->{_Loci}{$id}[$a]{'start'}<=>$self->{_Loci}{$id}[$b]{'start'}} $self->allKeys($id)){
		if ($pos>=$self->{_Loci}{$id}[$locus]{'start'} && $pos<=$self->{_Loci}{$id}[$locus]{'end'}){
			$index=$locus;
			last LOCUS;
			}
		elsif ($pos>$self->{_Loci}{$id}[$locus]{'end'}){next LOCUS}
		else {die "Could not find locus index for $id:$pos!\n"}
		}
	return $index;
	}


sub getreftoIndividuals{
	my ($self, $mode, $chrom, $id)=@_;
	if ($mode eq 'position'){ return($self->{_PosSummary}{$chrom}{$id}) }
	elsif ($mode eq 'locus'){ return($self->{_IndSummary}{$chrom}[$id]) }
	else {die "Mode not recognized!\n"}
	}

sub getreftoHaps{
	my ($self, $chrom, $id)=@_;
	return($self->{_HapSummary}{$chrom}[$id]);
	}


sub nextIterator{
	my $self=shift;
	return(each %{$self->{_Loci}});
	}


sub resetIterator{
	my $self=shift;
	keys %{$self->{_Loci}};
	return;
	}


sub addLocus{
	my ($self, %args)=@_;
	unless ( defined $args{'id'} && defined $args{'start'} && defined $args{'end'} ){ die "Locus location arguments missing!\n" }
	push @{ $self->{_Loci}{ $args{'id'} } }, { 'start'=>$args{'start'}, 'end'=>$args{'end'} };
	if ( defined $args{'refseq'} ){ $self->{_Loci}{ $args{'id'} }[-1]{'refseq'}=$args{'refseq'} }
	if ( defined $args{'info'} ){ $self->{_Loci}{ $args{'id'} }[-1]{'info'}=$args{'info'} }
	return;
	}


sub addRegions{
	my ($self, $ref_regions, $sort)=@_;
	my $lociadded=0;
	keys %{ $ref_regions->{_Loci} };
	CHROM: while ( my ($chrom, $ref_chrom)=each %{ $ref_regions->{_Loci} } ){
		if ( defined $self->{_Loci}{$chrom} ){ push @{ $self->{_Loci}{$chrom} }, @$ref_chrom }
		else { $self->{_Loci}{$chrom}=$ref_chrom }
		$lociadded+=scalar(@{ $ref_chrom });
		if ($sort){	# sort locus array by start position after adding new loci.
			my @sorted=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} };
			$self->{_Loci}{$chrom}=\@sorted;
			}
		}
	print STDERR "Added a total of $lociadded regions to the object.\n";
	return($lociadded);
	}


sub addRefseq{
	my ($self, $fasta)=@_;
	ID: while ( my ($id, $ref_id)=each %{ $self->{_Loci} } ){
		LOCUS: for my $ref_locus (@$ref_id){
			my $ref_seqstring=$fasta->randomAccess($id, $ref_locus->{'start'}, $ref_locus->{'end'}, 1);
			$ref_locus->{'refseq'}=$$ref_seqstring;
			print STDERR "Added reference sequence to $id:$ref_locus->{'start'}-$ref_locus->{'end'}.\n";
#			print STDERR "$ref_locus->{'refseq'}\n";
			}
		}
	return;
	}


sub readBED{
	my ($self, $filename, $format, $keepinfo)=@_;
	unless(@_==4){die "Wrong number of arguments!\n"}
	if ($filename eq "none"){ return }

	my %regions;

	open my $bed, "<", $filename or die "Could not open bed file $filename!\n";

	while (my $line= <$bed>){
		chomp $line;
		if ($line=~/^(\w+)\s+(\d+)\s+(\d+)\s*(.*)/){
			if ($2>$3){ warn "$1:$2-$3 has smaller end than start coordinate!\n"; next }
			if ($format eq 'bed'){push(@{$regions{$1}}, {'start'=>$2, 'end'=>($3-1)})}
			elsif ($format eq 'onebased'){push(@{$regions{$1}}, {'start'=>($2-1), 'end'=>($3-1)})}
			else {push(@{$regions{$1}}, {'start'=>$2, 'end'=>$3})}
			if ($keepinfo==1){$regions{$1}[-1]{'info'}=$4}
			}
		}
	close $bed;
	$self->{_Loci}=\%regions;
	return;
	}


sub readGTF{
	my ($self, $filename, $keepinfo, $type)=@_;
	unless(@_>=3){die "Wrong number of arguments!\n"}
	my %regions;

	open my $input, "<", $filename or die "Could not open bed file $filename!\n";

	while (my $line=<$input>){
		chomp $line;
		if ($line=~/^\w+\t\w+\t\w+/){
			my @col=split(/\t/, $line);
			unless (defined $type && $col[2] ne $type){
				my ($gene_id)=$col[8]=~/gene_id\s+"(\w+)";/;
				my ($transcript_id)=$col[8]=~/transcript_id\s+"(\w+)";/;
#				print STDERR "$col[0]:$col[3]-$col[4], $col[6], $gene_id\n";
				push( @{$regions{ $col[0] }}, {'type'=>$col[2], 'start'=>($col[3]-1), 'end'=>($col[4]-1), 'dir'=>$col[6], 'frame'=>$col[7], 'geneid'=>$gene_id, 'transcriptid'=>$transcript_id} );
				if ($keepinfo==1){ $regions{ $col[0] }[-1]{'info'}=$col[8] }
				}
			}
		}
	close $input;
	$self->{_Loci}=\%regions;
	return;
	}


sub readBins{
	my ($self, @filenames)=@_;
	unless(@_>=2){die "Wrong number of arguments!\n"}

	for my $filename (@filenames){
		open my $input, "<", $filename or die "Could not open bin file $filename!\n";

		while (my $line=<$input>){
			chomp $line;
			if ($line=~/^"?(\w+)"?\s+(\d+)/){
				if (exists $self->{_Bins}{$1}){
					if ($self->{_Bins}{$1} ne $2){ warn "GeneID $1 has already been assigned to a recombination bin! ", $self->{_Bins}{$1}, " vs. $2\n" }
					}
				else { $self->{_Bins}{$1}=$2 }
				}
			}
		close $input;
		}
	return;
	}


sub checkRegions{
	my ($self, $minlength)=@_;
	unless(@_==2){die "Wrong number of arguments!\n"}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $self->allKeys() ){
		my ($prevstart, $prevend);
		LOCUS: for my $locus ( sort { $a->{'start'}<=>$b->{'start'} } @{$self->{_Loci}{$chrom}} ){
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			my $locuslength=$end-$start+1;
			if ($locuslength<$minlength){ warn "Warning: locus $chrom:$start-$end has length below minimum locus length! $locuslength vs. $minlength\n" }
			if (defined $prevend && $start<=$prevend){ warn "Warning: overlap detected between $chrom:$prevstart-$prevend and $chrom:$start-$end!\n" }
			$prevstart=$start;
			$prevend=$end unless (defined $prevend && $prevend>$end);
			}
		}
	return;
	}

sub cleanRegions{
	my ($self, $minlength, $sort, $remove_overlap, $ref_bins)=@_;
	unless(@_>=3){die "Wrong number of arguments!\n"}
	my %cleaned;
	my %bins;
	%bins=map { $_=>1 } @$ref_bins if defined $ref_bins;

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $self->allKeys() ){
		my @loci;
		if ($sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{$self->{_Loci}{$chrom}} }
		else { @loci=@{$self->{_Loci}{$chrom}} }
		my $prevlocus;
		my $discard=0;
		LOCUS: for my $locus (@loci){
			$locus->{'bin'} //= $self->{_Bins}{ $locus->{'geneid'} };	# assign bin hash element if not already defined.
			if (defined $ref_bins && $ref_bins->[0] ne 'all'){
				next LOCUS unless ( defined $locus->{'bin'} && $bins{ $locus->{'bin'} } );
				}
			my $accepted=1;
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			my $locuslength=$end-$start+1;
			if ($locuslength<$minlength){
				warn "Warning: locus $chrom:$start-$end has length below minimum locus length! $locuslength vs. $minlength\n";
				next LOCUS;
				}
			unless (defined $prevlocus){ $prevlocus=$locus; next LOCUS }
			if ( $start<=$prevlocus->{'end'} ){
				my $infostring="$chrom:$prevlocus->{'start'}-$prevlocus->{'end'} (gene_id: $prevlocus->{'geneid'}, transcript_id: $prevlocus->{'transcriptid'}) and $chrom:$start-$end (gene_id: $locus->{'geneid'}, transcript_id: $locus->{'transcriptid'})";
#				print "Warning: overlap detected between $infostring!\n";
				if ( $prevlocus->{'geneid'} eq $locus->{'geneid'} ){
					unless ( $prevlocus->{'dir'} eq $locus->{'dir'} ){
						print STDERR "$infostring\n";
						die "Overlapping intervals have the same gene id but are on different strands!\n";
						}
					if ( $start==$prevlocus->{'start'} && $end==$prevlocus->{'end'} ){ next LOCUS }
					my $prevshift=$prevlocus->{'frame'};
					if ($prevlocus->{'dir'} eq '-'){ $prevshift=( $prevlocus->{'end'} - $prevlocus->{'start'} + 1 - $prevlocus->{'frame'} ) % 3 }
					my $shift=$locus->{'frame'};
					if ($locus->{'dir'} eq '-'){ $shift=( $locuslength - $locus->{'frame'} ) % 3 }
					if ( ($prevlocus->{'start'} + $prevshift - $start - $shift) % 3 ){	# check if both overlapping intervals are in the same reading frame.
						print STDERR "$infostring\n";
						warn "Overlapping intervals have the same gene id but are not in the same reading frame! Discarding both intervals.\n";
						print STDERR "$chrom:$prevlocus->{'start'}-$prevlocus->{'end'}, direction: $prevlocus->{'dir'}, frame: $prevlocus->{'frame'}, shift: $prevshift, codon start: ", $prevlocus->{'start'} + $prevshift, "\n"; 
						print STDERR "$chrom:$start-$end, direction: $locus->{'dir'}, frame: $locus->{'frame'}, shift: $shift, codon start: ", $start + $shift, "\n";
						$prevlocus={ %$prevlocus };	# make a copy of the prevlocus hash to avoid changing the original locus hash.
						$prevlocus->{'end'}=$end if $end > $prevlocus->{'end'};	# fuse intervals to detect further overlaps with either interval.
						$discard=1;
						next LOCUS;
						}
					else {
						print STDERR "$infostring\n";
						$prevlocus={ %$prevlocus };	# make a copy of the prevlocus hash to avoid changing the original locus hash.
						$prevlocus->{'end'}=$end if $end > $prevlocus->{'end'};
						warn "Overlapping intervals have same gene id, fusing intervals to $chrom:$prevlocus->{'start'}-$prevlocus->{'end'}!\n";
						next LOCUS;
						}
					}
				elsif ($remove_overlap){
					print STDERR "$infostring\n";
					warn "Overlapping regions have different gene ids, removing both intervals!\n";
					print STDERR "Removing $prevlocus->{'geneid'}\t$prevlocus->{'transcriptid'}\n";
					print STDERR "Removing $locus->{'geneid'}\t$locus->{'transcriptid'}\n";
					$prevlocus={ %$prevlocus };	# make a copy of the prevlocus hash to avoid changing the original locus hash.
					$prevlocus->{'end'}=$end if $end > $prevlocus->{'end'};	# fuse intervals to detect further overlaps with either interval.
					$discard=1;
					next LOCUS;
					}
				elsif ($prevlocus->{'dir'} ne $locus->{'dir'}){
					print STDERR "$infostring\n";
					warn "Overlapping regions have different gene ids and different directions, keeping both intervals!\n";
					}
				else {
					print STDERR "$infostring\n";
					warn "Overlapping regions have different gene ids and same directions, keeping longer interval!\n";
					if ( $locuslength < ($prevlocus->{'end'}-$prevlocus->{'start'}+1) ){ print "Removing $chrom:$start-$end!\n"; next LOCUS }
					else { print STDERR "Removing $chrom:$prevlocus->{'start'}-$prevlocus->{'end'}!\n"; $accepted=0 }
					}
				}
			$accepted=0 if $discard;
			push( @{$cleaned{$chrom}}, { %$prevlocus } ) if $accepted;
			$prevlocus=$locus;
			$discard=0;
			}
		push( @{$cleaned{$chrom}}, $prevlocus ) if (defined $prevlocus && !$discard);
		}
	my $cleaned_genes=Regions->new(\%cleaned);
	return($cleaned_genes);
	}

sub chopRegions{
	my ($self, $rmstart, $rmend, $minlength)=@_;
	unless(@_==4){die "Wrong number of arguments!\n"}
	my $minsize=$rmstart+$rmend+$minlength;

	CHROM: for my $chrom (keys %{$self->{_Loci}}){
		my @chopped;
		LOCUS: for my $locus ( sort { $a->{'start'}<=>$b->{'start'} } @{$self->{_Loci}{$chrom}} ){
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			my $locuslength=$end-$start+1;
			unless($locuslength){ warn "$chrom:$start-$end has smaller end than start coordinate!\n" }
			unless($locuslength>$minsize){ warn "$chrom:$start-$end is below minimum size after chopping, locus omitted!\n"; next LOCUS }
			$locus->{'start'}=$start+$rmstart;
			$locus->{'end'}=$end-$rmend;
			push @chopped, $locus;
			}
		$self->{_Loci}{$chrom}=\@chopped;
		}

	return;
	}


sub concatRegions{
	my ($self)=@_;
	unless(@_==1){die "Wrong number of arguments!\n"}

	my %concat;
	CHROM: for my $chrom (keys %{$self->{_Loci}}){
#		my $totsize=0;
		my ($ref_locus, $newstart, $prevend);
		LOCUS: for my $locus ( sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } ){
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			unless (defined $newstart){ $newstart=$start }
#			$totsize=$end-$newstart+1;
			if ( defined $prevend && $start>$prevend+1 ){
				push @{$concat{$chrom}}, { %$ref_locus };
				$concat{$chrom}[-1]{'start'}=$newstart;
				$concat{$chrom}[-1]{'end'}=$prevend;
				$newstart=$start;
#				$totsize=0;
				}
			$prevend=$end unless (defined $prevend && $prevend>$end);
			$ref_locus=$locus;
			}
		push @{$concat{$chrom}}, { %$ref_locus };
		$concat{$chrom}[-1]{'start'}=$newstart;
		$concat{$chrom}[-1]{'end'}=$prevend;
		}
	my $newloci=Regions->new(\%concat);
	return($newloci);
	}


sub windowingLoci{
	my ($self, $minlength, $windowsize, $stepsize, $offset, @chromosomes)=@_;
	my $verbose=0;
	unless(@_>=3){die "Wrong number of arguments!\n"}
	$stepsize=$windowsize unless (defined $stepsize);
	$offset=0 unless (defined $offset);
	unless (@chromosomes && $chromosomes[0] ne 'all'){ @chromosomes=keys %{$self->{_Loci}} }

	my %windows;
	CHROM: for my $chrom (@chromosomes){
		LOCUS: for my $locus (sort { $a->{'start'}<=>$b->{'start'} } @{$self->{_Loci}{$chrom}}){
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			my $locuslength=$end-$start+1;
			next LOCUS unless ($locuslength>=$minlength);
			my $windowstart=$start+$offset;
			print STDERR "$chrom, $start-$end\n" if $verbose;
			WINDOW: while ($windowstart<$end){
				my $windowend=$windowstart+$windowsize-1;
				if ($windowend>$end){
					$windowend=$end;
					last WINDOW if ($end-$windowstart+1<$minlength);	# Don't add incomplete windows if smaller than minlength.
					}
				push @{$windows{$chrom}}, {'start'=>$windowstart, 'end'=>$windowend};
#				print STDERR "$windowstart, $windowend\n";
				last WINDOW if $windowend==$end;
				$windowstart+=$stepsize;
				}
			}
		}
	my $newloci=Regions->new(\%windows);
	return($newloci);
	}


sub subsetLocibyID{
	my ($self, $ref_ids)=@_;
	unless(@_==2){die "Wrong number of arguments!\n"}
	my %subset;
	ID: for my $id (@$ref_ids){
		unless (defined $self->{_Loci}{$id} ){ warn "Could not find $id in regions object!\n"; next ID }
		push @{ $subset{$id} }, { %$_ } foreach (@{ $self->{_Loci}{$id} });
		}
	my $newloci=Regions->new(\%subset);
	return($newloci);
	}

sub subsetLoci{	# Attention: Uses hash references. Any change to the subsetted loci will also affect the full set of loci!
	my ($self, $maxsize, @chromosomes)=@_;
	unless(@_>=2){die "Wrong number of arguments!\n"}
	unless (@chromosomes && $chromosomes[0] ne 'all'){ @chromosomes=sort {Misc::expand($a) cmp Misc::expand($b)} keys %{ $self->{_Loci} } }
	my %subset;
	my $prevend=0;
	my $rt=0;
	CHROM: for my $chrom (@chromosomes){
		next CHROM unless (defined $self->{_Loci}{$chrom} );
		unless ( defined $self->{_LastIndex}{$chrom} ){
			@{$self->{_Loci}{$chrom} }=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} };	# make sure to have loci sorted by scaffold position.
			$self->{_LastIndex}{$chrom}=0;
			}
		my $endindex=$self->getSize($chrom)-1;
		next CHROM if $self->{_LastIndex}{$chrom}>$endindex;
		LOCUS: for my $locus ( $self->{_LastIndex}{$chrom}..$endindex){
#			if ($self->getExcluded($chrom, $locus)){ next LOCUS }
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
			if ($start<$prevend){ $rt+=$end-$prevend }
			else { $rt+=$end-$start+1 }
			if ($rt<=$maxsize){
				push @{$subset{$chrom}}, $self->{_Loci}{$chrom}[$locus];	# makes sure to take over all the information in the hash.
#				$self->setExcluded($chrom, $locus, 1);
				}
			else {$self->{_LastIndex}{$chrom}=$locus; last CHROM}
			$prevend=$end;
			}
		$self->{_LastIndex}{$chrom}=$self->getSize($chrom);
		last CHROM if exists $subset{$chrom};	# only process one scaffold at a time to prevent disruption of scaffold order on chromosome (scaffolds are stored as hash elements)!
		}
	my $newloci=Regions->new(\%subset);
	return($newloci);
	}


sub excludeRegions{
	my ($self, $regions, $min_dist, $verbose)=@_;

	my %reduced;	

	SCAFF: for my $chrom (sort keys %{$self->{_Loci}}){
		LOCUS: for my $locus ( @{ $self->{_Loci}{$chrom} } ){
			if ( defined $regions->getValue($chrom) ){
				REG: for my $reg ( sort { $a->{'start'}<=>$b->{'start'} } $regions->allValues($chrom) ){
					if ( ($reg->{'end'}+$min_dist)>=$locus->{'start'} && ($reg->{'start'}-$min_dist)<=$locus->{'end'} ){
						warn "Overlap detected between $chrom $reg->{'start'}-$reg->{'end'} and $locus->{'start'}-$locus->{'end'}\n" if $verbose;
						next LOCUS;
						}
					elsif ( $reg->{'start'}>$locus->{'end'} ){ last REG }
					}
				}
			push @{ $reduced{$chrom} }, $locus;
			}
		}

	my $newloci=Regions->new(\%reduced);
	return($newloci);
	}


sub splitRegions{  # needs more testing!
	my ($self, $regions, $min_dist, $verbose)=@_;
	my %splitted;	

	SCAFF: for my $chrom (sort keys %{$self->{_Loci}}){
		LOCUS: for my $locus ( @{ $self->{_Loci}{$chrom} } ){
			if ( defined $regions->getValue($chrom) ){
				REG: for my $reg (sort { $a->{'start'}<=>$b->{'start'} } $regions->allValues($chrom)){
				    my $regstart=$reg->{'start'}-$min_dist;
				    my $regend=$reg->{'end'}+$min_dist;
				    if ($regstart<=$locus->{'start'}){
				        if ($regend<$locus->{'start'}){ next REG }
				        elsif ($regend>=$locus->{'end'}){ next LOCUS } # locus is completely nested within exlusion region.
				        elsif ($regend<$locus->{'end'}){ $locus->{'start'}=$regend+1 } # left-side overlap.
				        }
				    else {
				        if ($regend<$locus->{'end'}){    # region is completely nested within locus, need to split.
                            warn "Split required for locus $locus->{'start'}-$locus->{'end'} due to excluded region $chrom $reg->{'start'}-$reg->{'end'}.\n" if $verbose;
				            push @{ $splitted{$chrom} }, { %$locus };
				            $splitted{$chrom}[-1]{'end'}=$regstart-1;
				            $locus->{'start'}=$regend+1;
				            }
				        elsif ($regstart<=$locus->{'end'}){ $locus->{'end'}=$regstart-1 } # right-side overlap.
				        else { last REG }
				        }
					}
				}
			push @{ $splitted{$chrom} }, $locus;
			}
		}

	my $newloci=Regions->new(\%splitted);
	return($newloci);
	}


sub testProximity{
	my ($self, $candidates, $min_dist, $verbose)=@_;

	my @accepted;

	CAND: for my $cand ($candidates->allValues()){
		my $id=$cand->{'id'};
		my $start=$cand->{'start'};
		my $end=$cand->{'end'};

		if (defined $self->{_Loci}{$id}){
			REG: for my $reg (@{$self->{_Loci}{$id}}){
				if (($reg->{'end'}+$min_dist)>=$start && ($reg->{'start'}-$min_dist)<=$end){
					warn "Overlap detected between $id: $reg->{'start'}-$reg->{'end'} and $start-$end, locus rejected!\n" if $verbose;
					push @accepted, 0;
					next CAND;
					}
				}
			}
		print STDERR "No overlap detected for locus $id: $start-$end, locus accepted.\n" if $verbose;
		
		$self->addLocus('id'=>$id, 'start'=>$start, 'end'=>$end);	# uses now named arguments.
		push @accepted, 1;
		}

	return;
	}


sub generateFasta{
	my ($self, $ref_fasta, $ref_fastaarray, $ref_genotypes, $ref_missing, $ref_minind, $groups, $ref_mingroupind, $mincov, $excluderefN, $haplotize, $ref_subind)=@_;

	my $verbose=0;
	my @individuals;
	my (@popsizes, @popsum, @groupsizes, @groupsum);
	my $npops=$ref_missing->getNPop();
	my $totsize=0;
	for my $group ( 0..$#{$groups} ){
		push @groupsum, $totsize;
		my $groupsize=0;
		for my $pop ( @{$groups->[$group]} ){
			push @popsum, $totsize;
			my $temp=$ref_missing->getNInd($pop);
			push @popsizes, $temp;
			if (defined $ref_subind){ push @individuals, ( map { $totsize+$_ } @{ $ref_subind->[$pop] } ) }
			else { push @individuals, ($totsize..$totsize+$temp-1) }
			$totsize+=$temp;
			$groupsize+=$temp;
			print STDERR "pop $pop: popsize $temp, popsum $popsum[$pop], totsize: $totsize\n" if $verbose;
			}
		push @groupsizes, $groupsize;
		print STDERR "group $group: groupsize $groupsize, groupsum $groupsum[$group], totsize: $totsize\n" if $verbose;
		}
	die "Wrong length of array of fasta objects provided!\n" unless (@$ref_fastaarray==@individuals);
	print STDERR join(',', @individuals), "\n" if $verbose;

	CHROM: for my $chrom ( $self->allKeys() ){
		my $prevend;
		LOCUS: for my $locus ( sort { $self->{_Loci}{$chrom}[$a]{'start'}<=>$self->{_Loci}{$chrom}[$b]{'start'} } $self->allKeys($chrom) ){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
			if (defined $prevend){
				if ($start-$prevend>1){ $start=$prevend+1; warn "Gap detected between two consecutive loci on the same chromosome. Start changed to $start.\n" }
				elsif (!$start-$prevend){ $start=$prevend+1; warn "Overlap detected between two consecutive loci on the same chromosome. Start changed to $start.\n"}
				}
			my $locuslength=$end-$start+1;
			my @seq;
			my $ref_seqstring;
			if (defined $ref_fasta){ $ref_seqstring=$ref_fasta->randomAccess($chrom, $start, $end, 1) }
			else { my $seqstring="A" x $locuslength; $ref_seqstring=\$seqstring }
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length at $chrom:$start-$end. ", length($$ref_seqstring), " vs. $locuslength!\n" }
#			for my $allind (0..$totsize-1){ $seq[$allind]=$$ref_seqstring }	# set the sequences of all individuals to reference.
			my ($totalpos, $hardmasked)=(0, 0);
			my $basestring;
			POSITION: for my $pos ($start..$end){
				$basestring="";
				my $refbase=substr($$ref_seqstring, 0, 1, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				my $offset=$pos-$start;
				if ($excluderefN && $refbase eq 'N'){ ++$hardmasked; $basestring="N" x $totsize; next POSITION }	# discard if base in reference sequence is hardmasked.
				my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
				my $groupbases;
				GROUP: for my $group ( 0..$#{$groups} ){
					$groupbases="";
					my $validgroupind=0;
					POP: for my $pop ( @{$groups->[$group]} ){
						my $validind=0;
						IND: for my $ind (0..$popsizes[$pop]-1){
							my $allind=$ind+$popsum[$pop];
							my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
							unless (defined $coverage && $coverage>=$mincov){ $groupbases.="N"; next IND }	# missing data, replace reference base at offset with N.
							if ($snp){	# if position not variable but covered, retain reference base.
								my $index=2*$ind;
								my $allind=$ind+$popsum[$pop];
								my $allele1=$ref_genotypes->getBase($chrom, $pos, $pop, $index);
								unless (defined $allele1 && $allele1 ne "N"){ $groupbases.="N"; next IND }
								my $allele2=$ref_genotypes->getBase($chrom, $pos, $pop, ++$index);
								unless (defined $allele2 && $allele2 ne "N"){ $groupbases.="N"; next IND }
								if ($allele1 eq $allele2){ $groupbases.=$allele1 }
								elsif ($haplotize){ $groupbases.=($allele1, $allele2)[ int(rand(2)) ] } 
								else { $groupbases.=Misc::iupac($allele1, $allele2) }
								}
							else { $groupbases.=$refbase }
							++$validind
							}
						if ($validind<$ref_minind->[$pop]){	# not enough individuals callable in a population, set bases to N for all individuals in the group;
							$groupbases="N" x $groupsizes[$group];
							next GROUP;
							}
						$validgroupind+=$validind;
						}
					if ($validgroupind<$ref_mingroupind->[$group]){	# not enough individuals callable in a group, set bases to N for all individuals in the group;
						$groupbases="N" x $groupsizes[$group];
						next GROUP;
						}
					die "Wrong length of groupbases string at $chrom:$pos, ", length($groupbases), " vs. $groupsizes[$group]!\n$groupbases\n" unless ( length($groupbases)==$groupsizes[$group] );
					} continue { $basestring.=$groupbases }
				die "Wrong length of basestring at $chrom:$pos, ", length($basestring), " vs. $totsize!\n$basestring\n" unless (length($basestring)==$totsize);
				} continue { $seq[$_].=substr($basestring, 0, 1, "") foreach (0..$totsize-1) }
			for my $selindex (0..$#individuals){
				my $selind=$individuals[$selindex];
				unless ( length($seq[$selind])==$locuslength){ die "Wrong length of sequence string $chrom:$start-$end of individual $selind, ",
					length($seq[$selind]), " vs. $locuslength, concatenation aborted!\n" }
				$ref_fastaarray->[$selindex]->concatSeq($chrom, \$seq[$selind]);
				}
			$prevend=$end;
			}
		}
	return;
	}


sub printLoci{
	my ($self, $outfile, $printinfo)=@_;
	my $output;

	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", $outfile or die "Could not open output file!\n"}

	SCAFF: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		LOCUS: for my $locus (sort {$a->{'start'}<=>$b->{'start'}} @{$self->{_Loci}{$chrom}}){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1;
			print $output "\t", $locus->{'info'} || '' if $printinfo;
			print $output "\n";
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}

	return;
	}

sub printLociAttributes{
	my ($self, $outfile, $ref_attributes, $sort, $append)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT }
	elsif ($append){ open $output, ">>", $outfile or die "Could not open output file!\n" }
	else { open $output, ">", $outfile or die "Could not open output file!\n" }

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if ($sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1, "\t", $locus->{'info'} || "";
			if (ref($locus->{$ref_attributes->[0]} ) eq "ARRAY"){
				for my $group (0..$#{ $locus->{$ref_attributes->[0]} }){
					print $output "\t$locus->{$_}[$group]" foreach (@$ref_attributes);
					}
				}
			else { print $output "\t$locus->{$_}" foreach (@$ref_attributes); }
			print $output "\n";
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return;
	}


1;

