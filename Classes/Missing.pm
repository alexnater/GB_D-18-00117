package Missing;

use strict;
use warnings;
use Classes::Genotypes;
use Classes::Misc;


sub new{
	my ($class, $allsiteslist, $type, $ref_sampleids)=@_;
	my $self={};

	$self->{_Filetype}=defined $type ? $type : 'vcf';
	my $coloffset=(defined $type && $type eq 'depth') ? 1 : 8;

	if (defined $allsiteslist){
		my ($ref_indnames, $ref_files, $ref_samples)=Misc::getIndnames($allsiteslist, $coloffset);
        ($ref_samples, $ref_indnames)=Misc::getIndicesfromIndnames($ref_indnames, $ref_samples, $ref_sampleids) if (defined $ref_sampleids);
		$self->{_Files}=$ref_files;
		$self->{_Samples}=$ref_samples;
		$self->{_Indnames}=$ref_indnames;
       	$self->{_Popsizes}=[ map { scalar(@$_) } @$ref_samples ];
	  	$self->{_Popsums}=Misc::getPopsums($self->{_Popsizes});
       	$self->{_Popindices}=Misc::getPopindices($ref_indnames);
        $self->{_Allindmap}=Misc::getAllindMap($self->{_Popsizes}, $self->{_Popindices});
       	$self->{_Simulated}=0;
		}
    else {
        $self->{_Simulated}=1;
        }

	bless $self, $class;
	return $self;
	}


sub setPopulations{
	my ($self, $ref_populations, $verbose)=@_;
	unless (@_>=2){ die "Wrong number of arguments!\n" }

	# get the sample indices from the vcf files:
	my $ref_sampleids=$ref_populations->getSampleids();
    my ($ref_subsamples, $ref_subindnames)=Misc::getIndicesfromIndnames($self->{_Indnames}, $self->{_Samples}, $ref_sampleids);
    
    # check if all samples in the poplist have been found in the vcf files:
    my @vcf_sampleids=map { @$_ } @$ref_subindnames;
    if (@vcf_sampleids!=@$ref_sampleids){
        die "ERROR: Not all samples in the poplist have been found in the missing data files!\n";
#        $ref_populations->setSubset(['all'], \@vcf_sampleids);
        }

    if ($verbose){
        print STDERR "File\tIndividualID\tGrouplabel\tGroupindex\tPoplabel\tPopindex\n";
        for my $file (0..$#{ $ref_subindnames }){
            for my $ind (0..$#{ $ref_subindnames->[$file] }){
                my $indname=$ref_subindnames->[$file][$ind];
                my ($poplabel, $popidx)=$ref_populations->getPop($indname);
                my ($grouplabel, $groupidx)=$ref_populations->getGroup($poplabel);
                print STDERR "$file\t$indname\t$grouplabel\t$groupidx\t$poplabel\t$popidx\n";
                }
            }
        }
    
    $self->{_Samples}=$ref_subsamples;
    $self->{_Indnames}=$ref_populations->getSampleids_bypop();
  	$self->{_Popsizes}=$ref_populations->getPopsizes();
  	$self->{_Popsums}=Misc::getPopsums($self->{_Popsizes});
  	$self->{_Groups}=$ref_populations->getGroups();
   	$self->{_Popindices}=Misc::getPopindices($ref_subindnames, $ref_populations->getPoplabels(), $ref_populations->getPopmap());
	$self->{_Popfiles}=Misc::getPopfiles($self->{_Popindices});
    $self->{_Allindmap}=Misc::getAllindMap($self->{_Popsizes}, $self->{_Popindices});
    return;
    }


sub getValuebyPop{
	my ($self, $id, $startpos, $pop, $ind, $offset)=@_;

	if ( $self->{_Simulated} ){ return 255 }
	elsif (@_==5){ return($self->{_Missing}{$id}{$startpos}[$self->{_Popsums}[$pop]+$ind]) }
	elsif (@_==6){
		my $substring=substr($self->{_Missing}{$id}{$startpos}[$self->{_Popsums}[$pop]+$ind], $offset, 1);
		unless (defined $substring){ warn "Invalid access: $id, $startpos, $pop, $ind, $offset!\n" }
		return( ord($substring) );
		}
	else { die "Invalid number of arguments!\n" }
	}

sub getValue{
	my ($self, $id, $startpos, $allind, $offset)=@_;

	if ( $self->{_Simulated} ){ return 255 }
	elsif (@_==4){ return($self->{_Missing}{$id}{$startpos}[$allind]) }
	elsif (@_==5){
		my $substring=substr($self->{_Missing}{$id}{$startpos}[$allind], $offset, 1);
		unless (defined $substring){ warn "Invalid access: $id, $startpos, $allind, $offset!\n" }
		return( ord($substring) );
		}
	else { die "Invalid number of arguments!\n" }
	}


sub setValue{
	my ($self, $id, $startpos, $allind, $offset, $value)=@_;
	unless (@_==6){die "Wrong number of arguments!\n"}

	if ($value>255){$value=255}
	substr($self->{_Missing}{$id}{$startpos}[$allind], $offset, 1, chr($value) );
	return($value);
	}


sub getSize{
	my ($self, $id, $pos, $pop)=@_;

	if (@_==1){ return(scalar(keys %{$self->{_Missing}})) }
	elsif (@_==2){ return(scalar(keys %{$self->{_Missing}{$id}})) }
	elsif (@_==3){ return(scalar @{$self->{_Missing}{$id}{$pos}}) }
	elsif (@_==4){ return(scalar @{$self->{_Missing}{$id}{$pos}[$pop]}) }
	else { die "Invalid number of arguments!\n" }
	}


sub getPosition{
	my ($self, $id)=@_;
	unless (@_==2){die "Wrong number of arguments!\n"}

	my $position=each %{$self->{_Missing}{$id}};
	keys %{$self->{_Missing}{$id}};

	return($position);
	}


sub getGroups{
	my ($self, $group)=@_;
	if (@_==1){ return($self->{_Groups}) }
	elsif (@_==2){ return($self->{_Groups}[$group]) }
	else { die "Wrong number of arguments!\n" }
	}

sub getGroupsizes{
	my ($self, $group)=@_;
	if (@_==1){
		return([ map { Misc::sum(@{ $self->{_Popsizes} }[@$_]) } @{ $self->{_Groups} } ]);
		}
	elsif (@_==2){
		return(Misc::sum(@{ $self->{_Popsizes} }[@{ $self->{_Groups}[$group] }]) );
		}
	else { die "Wrong number of arguments!\n" }
	}

sub getNPop{
	unless (@_==1){die "Wrong number of arguments!\n"}
	my $self=shift;
	return(scalar @{$self->{_Popsizes}});
	}

sub getNInd{
	my ($self, $pop)=@_;
	if (@_==1){ return($self->{_Popsizes}) }
	elsif (@_==2){ return($self->{_Popsizes}[$pop]) }
	else { die "Wrong number of arguments!\n" }
	}

sub getPopsums{
	unless (@_==1){die "Wrong number of arguments!\n"}
	my $self=shift;
	return($self->{_Popsums});
	}

sub getIndnames{
	my ($self, $pop, $ind)=@_;
	if (@_==1){ return($self->{_Indnames}) }
	elsif (@_==3){ return($self->{_Indnames}[$pop][$ind]) }
	else { die "Wrong number of arguments!\n" }
	}


sub getMQ{
	my ($self, $id, $startpos, $pop, $offset)=@_;
	if ( $self->{_Simulated} ){ return 255 }

	my $file=$self->{_Popfiles}{$pop};
	if (@_==4){ return($self->{_MQ}{$id}{$startpos}[$file]) }
	elsif (@_==5){
		my $substring=substr($self->{_MQ}{$id}{$startpos}[$file], $offset, 1);
		unless (defined $substring){ warn "Invalid access: $id, $startpos, $pop, $offset!\n" }
		return( ord($substring) );
		}
	else { die "Invalid number of arguments!\n" }
	}


sub readSamtoolsDepth{	
	my ($self, $ref_loci, $ref_popsizes, $excludeRefN, $verbose)=@_;

    $self->{_Missing}=undef;    # reset storage of missing data.
	unless (defined $self->{_Files} && defined $self->{_Samples} ){ die "List and folder of depth files have not been defined!\n" }
	my @popsizes=@{ $self->{_Popsizes} };
    my @filesizes=map { scalar(@$_) } (@{ $self->{_Samples} });

	if (defined $ref_popsizes){	# allows for population size check.
		if (@$ref_popsizes!=@popsizes){ warn "Wrong number of populations specified (", scalar(@$ref_popsizes), " vs. ", scalar(@popsizes), ")!\n" }
		for my $pop (0..$#{ $ref_popsizes } ){
			if ( $ref_popsizes->[$pop]!=$popsizes[$pop] ){ warn "Wrong number of samples specified for population $pop ($ref_popsizes->[$pop] vs. ", scalar($popsizes[$pop]), ")!\n" }
			}
		}

	CHROM: for my $chrom ( $ref_loci->allKeys() ){
		LOCUS: for my $locus ( $ref_loci->allKeys($chrom) ){
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			my $locuslength=$end-$start+1;
			my $reference=$ref_loci->getValue($chrom, $locus, 'refseq');
			unless (defined $reference){ warn "Reference sequence for $chrom:$start-$end is not defined!\n" if $verbose; $reference="" }
			my @missing;

			FILE: for my $file (0..$#{ $self->{_Files} }){
			    my @samples=@{ $self->{_Samples}[$file] };
				my @allindmap=@{ $self->{_Allindmap}[$file] };
                unless (@allindmap==@samples){ die "Wrong length of population assignment map for depth file $self->{_Files}[$file]!\n$!\n" }
				unless (@samples){
				    warn "No samples to process for file $file!\n";
				    next FILE;    # skip file if no samples to process.
				    }

				open my $vcffile, "-|", "tabix $self->{_Files}[$file] $chrom:$start-$end" or die "Could not open depth file $self->{_Files}[$file]!\n$!\n";
				print STDERR "Acquiring coverage data for $chrom:$start-$end, depth file $file.\n" if $verbose;
				my $refseq=$reference;	# make hard copy of reference sequence.
				my $prevpos=$start-1;
				my ($scaffold, $pos, $mq, $vq, $dp);
				LINE: while (<$vcffile>){
					if(/^[\w\-:]+\s/){
						my @line=split("\t", $_);
						$scaffold=$line[0];
						$pos=$line[1];

						if ($pos-$prevpos!=1){
							my $gapsize=$pos-$prevpos-1;
							for my $allind (@allindmap){ $missing[$allind].=(chr(0) x $gapsize) }
							substr($refseq, 0, $gapsize, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
							}
						$prevpos=$pos;

						my $refbase=uc( substr($refseq, 0, 1, "") );
						if ($excludeRefN){
							if ($refbase eq 'N' ){	# discard if base in reference sequence is hardmasked.
								for my $allind (@allindmap){ $missing[$allind].=chr(0) }
								next LINE;
								}
							}

						for my $idx (0..$#samples){
						    my $allind=$allindmap[$idx];
							my $col=$samples[$idx]+1;
							if ($col>$#line){ die "Invalid sample index specified for depth file $self->{_Files}[$file]!\n" }
							my $coverage=int($line[$col]);
							if (defined $coverage && $coverage<255){ $missing[$allind].=chr($coverage) }
							elsif (defined $coverage){ $missing[$allind].=chr(255) }
							else { $missing[$allind].=chr(0) }
							}
						}
					}
				close $vcffile;
				my $gapsize=defined $pos ? $end-$pos: $locuslength;
				for my $allind (@allindmap){
					if ($gapsize){ $missing[$allind].=(chr(0) x $gapsize) }
					my $stringlength=length( $missing[$allind] );
					if ($stringlength != $locuslength){
						warn "Non-matching string length for $chrom:$start-$end! $stringlength vs. $locuslength\n";
						}
					}
				}
			$self->{_Missing}{$chrom}{$start-1}=\@missing;			
			}
		}
	return;
	}


sub readAllsites{	
	my ($self, $ref_loci, $minmq, $minvq, $ref_popsizes, $excludeRefN, $recordmq, $verbose)=@_;

    # reset data:
	$self->{_Missing}=undef;
	$self->{_MQ}=undef;
	
	unless (defined $self->{_Files} && defined $self->{_Samples} ){ die "List and folder of depth files have not been defined!\n" }
	my @popsizes=@{ $self->{_Popsizes} };
    my @filesizes=map { scalar(@$_) } (@{ $self->{_Samples} });

	if (defined $ref_popsizes){	# allows for population size check.
		if (@$ref_popsizes!=@popsizes){ warn "Wrong number of populations specified (", scalar(@$ref_popsizes), " vs. ", scalar(@popsizes), ")!\n" }
		for my $pop (0..$#{ $ref_popsizes } ){
			if ( $ref_popsizes->[$pop]!=$popsizes[$pop] ){ warn "Wrong number of samples specified for population $pop ($ref_popsizes->[$pop] vs. ", scalar($popsizes[$pop]), ")!\n" }
			}
		}

	CHROM: for my $chrom ( $ref_loci->allKeys() ){
		LOCUS: for my $locus ( $ref_loci->allKeys($chrom) ){
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			my $locuslength=$end-$start+1;
			my $reference=$ref_loci->getValue($chrom, $locus, 'refseq');
			unless (defined $reference){ warn "Reference sequence for $chrom:$start-$end is not defined!\n" if $verbose; $reference="" }
			my (@missing, @mq_by_file);

			FILE: for my $file (0..$#{ $self->{_Files} }){
			    my @samples=@{ $self->{_Samples}[$file] };
				my @allindmap=@{ $self->{_Allindmap}[$file] };
                unless (@allindmap==@samples){ die "Wrong length of population assignment map for allsites file $self->{_Files}[$file]!\n$!\n" }
				unless (@samples){
				    warn "No samples to process for file $file!\n";
				    next FILE;    # skip file if no samples to process.
				    }

				open my $vcffile, "-|", "tabix $self->{_Files}[$file] $chrom:$start-$end" or die "Could not open allsites vcf file $self->{_Files}[$file]!\n$!\n";
				print STDERR "Acquiring coverage data for $chrom:$start-$end, allsites file $file.\n" if $verbose;
				my $refseq=$reference;	# make hard copy of reference sequence.
				my $prevpos=$start-1;
				my ($scaffold, $pos, $mq, $vq, $dp);
				LINE: while (<$vcffile>){
					if(/^[\w\-:]+\s/){
						my @line=split("\t", $_);
						$scaffold=$line[0];
						$pos=$line[1];

						if ($pos-$prevpos!=1){
							my $gapsize=$pos-$prevpos-1;
							for my $allind (@allindmap){ $missing[$allind].=(chr(0) x $gapsize) }
							$mq_by_file[$file].=(chr(0) x $gapsize) if ($recordmq);
							substr($refseq, 0, $gapsize, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
							}
						$prevpos=$pos;

						my $refbase=uc( substr($refseq, 0, 1, "") );
						if ($excludeRefN){
							if ($refbase eq 'N' ){	# discard if base in reference sequence is hardmasked.
								for my $allind (@allindmap){ $missing[$allind].=chr(0) }
    							$mq_by_file[$file].=chr(0) if ($recordmq);
								next LINE;
								}
							}

						($mq)=($line[7]=~/MQ=(\d+\.?\d*)/);
						if ($recordmq){
							if (defined $mq && $mq<255){ $mq_by_file[$file].=chr( int($mq) ) }
							elsif (defined $mq){ $mq_by_file[$file].=chr(255) }
							else { $mq_by_file[$file].=chr(0) }
							}

						$vq=$line[5];
						unless ($vq ne '.' && $vq>=$minvq){
							for my $allind (@allindmap){ $missing[$allind].=chr(0) }
							next LINE;
							}
						unless (defined $mq && $mq>=$minmq){
							for my $allind (@allindmap){ $missing[$allind].=chr(0) }
							next LINE
							}

						$dp=$line[4] eq '.' ? 1 : 2;

						for my $idx (0..$#samples){
						    my $allind=$allindmap[$idx];
							my $col=$samples[$idx]+8;
							if ($col>$#line){ die "Invalid sample index specified for depth file $self->{_Files}[$file]!\n" }
							my @data=split(":", $line[$col]);
							if (@data>1 && $data[$dp]<255){ $missing[$allind].=chr($data[$dp]) }
							elsif (@data>1){ $missing[$allind].=chr(255) }
							else { $missing[$allind].=chr(0) }
							}
						}
					}
				close $vcffile;
				my $gapsize=defined $pos ? $end-$pos: $locuslength;
				if ($recordmq && $gapsize){ $mq_by_file[$file].=(chr(0) x $gapsize) }
				for my $allind (@allindmap){
					if ($gapsize){ $missing[$allind].=(chr(0) x $gapsize) }
					my $stringlength=length($missing[$allind]);
					if ($stringlength != $locuslength){
						warn "Non-matching string length for $chrom:$start-$end! $stringlength vs. $locuslength\n";
						}
					}
				}
			$self->{_Missing}{$chrom}{$start-1}=\@missing;
			$self->{_MQ}{$chrom}{$start-1}=\@mq_by_file;			
			}
		}
	return;
	}


1;

