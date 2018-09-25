package Genotypes;

use strict;
use warnings;
use Classes::Misc;
use File::Basename;


sub new{
	my ($class, $vcflist, $ref_sampleids)=@_;
	my $self={};

	if (defined $vcflist){
		my ($ref_indnames, $ref_files, $ref_samples)=Misc::getIndnames($vcflist);
        ($ref_samples, $ref_indnames)=Misc::getIndicesfromIndnames($ref_indnames, $ref_samples, $ref_sampleids) if (defined $ref_sampleids);
		$self->{_Files}=$ref_files;
		$self->{_Samples}=$ref_samples;
		$self->{_Indnames}=$ref_indnames;
       	$self->{_Popsizes}=[ map { scalar(@$_) } @$ref_samples ];
       	$self->{_Popindices}=Misc::getPopindices($ref_indnames);
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
        die "ERROR: Not all samples in the poplist have been found in the genotype files!\n";
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
  	$self->{_Groups}=$ref_populations->getGroups();
   	$self->{_Popindices}=Misc::getPopindices($ref_subindnames, $ref_populations->getPoplabels(), $ref_populations->getPopmap());
    return;
    }


sub getValue{
	my ($self, $chrom, $pos, $pop, $index)=@_;
	if (@_==3){ return($self->{_Genotypes}{$chrom}{$pos}) }
	elsif (@_==4){ return($self->{_Genotypes}{$chrom}{$pos}[$pop]) }
	elsif (@_==5){
		return( substr($self->{_Genotypes}{$chrom}{$pos}[$pop], $index, 1) );
		}
	else { die "Invalid number of arguments!\n" }
	}

sub setValue{
	unless (@_==6){die "Wrong number of arguments!\n"}
	my ($self, $chrom, $pos, $pop, $index, $value)=@_;
	substr($self->{_Genotypes}{$chrom}{$pos}[$pop], $index, 1, $value);
	return($value);
	}


sub allKeys{
	unless (@_==2){die "Wrong number of arguments!\n"}
	my ($self, $id)=@_;
	return( keys %{$self->{_Genotypes}{$id}} )
	}


sub getAncestral{
	unless (@_==3){die "Wrong number of arguments!\n"}
	my ($self, $id, $pos)=@_;
	if (exists $self->{_Ancestral}{$id}{$pos}){ return($self->{_Ancestral}{$id}{$pos}) }
    else { return(4) }
	}

sub setAncestral{
	unless (@_==4){ die "Wrong number of arguments!\n" }
	my ($self, $id, $pos, $value)=@_;
	$self->{_Ancestral}{$id}{$pos}=$value;
	return;
	}


sub getBase{
	my ($self, $chrom, $pos, $i, $j)=@_;
	if (@_==3){ return( split('', $self->{_Bases}{$chrom}{$pos}) ) }	# returns all bases present at site as list.
	elsif (@_==4){
		if ($i eq '.'){ return("N") } 
		elsif ($i>=0 && $i<4){ return( substr($self->{_Bases}{$chrom}{$pos}, $i, 1) ) }
		else { return }
		}
	elsif (@_==5){
		my $gt=substr($self->{_Genotypes}{$chrom}{$pos}[$i], $j, 1);
		unless (defined $gt){ return }
		if ($gt eq '.'){ return("N") }
		else { return( substr($self->{_Bases}{$chrom}{$pos}, $gt, 1) ) }
		}
	else {die "Invalid number of arguments!\n"}
	}

sub getAlleleforBase{
	unless (@_==3 || @_==4){ die "Wrong number of arguments!\n" }
	my ($self, $chrom, $pos, $base)=@_;
	my @baselist=split('', $self->{_Bases}{$chrom}{$pos});
	my %basehash=map { $baselist[$_] => $_ } (0..$#baselist);
	if (defined $base){ return( $basehash{$base} ) }
	else { return(\%basehash) }		
	}


sub isAllele{
	unless (@_==4){die "Wrong number of arguments!\n"}
	my ($self, $id, $pos, $allele)=@_;
	unless ($allele>=0 && $allele<=3){die "Invalid allele specifier!\n"}
	if ($self->{_Alleles}{$id}{$pos} & (1 << $allele)){return 1}	# shift first bit to the left to position specified by $allele and test if bit is set.
	else {return 0}
	}

sub getAlleles{
	unless (@_==3){die "Wrong number of arguments!\n"}
	my ($self, $id, $pos)=@_;
	my @alleles=map { $self->{_Alleles}{$id}{$pos} & (1 << $_) ? $_ : () } (0..3);
	return (@alleles);
	}

sub getPolarizedAlleles{	# returns a list where the first element is the ancestral allele, even if not physically present in the set.
	unless (@_==3 || @_==4){die "Wrong number of arguments!\n"}
	my ($self, $id, $pos, $ancestral)=@_;
	$ancestral=$self->{_Ancestral}{$id}{$pos} unless (defined $ancestral);
	my $bits=$self->{_Alleles}{$id}{$pos};	# make copy of bitstring.
	if (defined $ancestral && $ancestral>=0 && $ancestral<4){
		$bits &= ~(1 << $ancestral);	# unset bit for ancestral allele.
		}
	else { $ancestral=undef }
	my @alleles=map { $bits & (1 << $_) ? $_ : () } (0..3);
	return ($ancestral, @alleles);
	}

sub setAlleles{
	unless (@_>=4){die "Wrong number of arguments!\n"}
	my ($self, $id, $pos, @alleles)=@_;
	for my $allele (@alleles){
		unless ($allele>=0 && $allele<=3){die "Invalid allele specifier!\n"}
		$self->{_Alleles}{$id}{$pos} |= (1 << $allele);	# shift first bit to the left to position specified by $allele and set bit to 1.
		}
	return;
	}

sub NAlleles{
	unless (@_==3){die "Wrong number of arguments!\n"}
	my ($self, $id, $pos)=@_;
	my $count1=unpack( '%32b*', chr($self->{_Alleles}{$id}{$pos}) );
	return($count1);
	}


sub getCoverage{
	my ($self, $id, $pos, $pop, $ind)=@_;
	if ( $self->{_Simulated} ){ return 255 }
	elsif (@_==3){ return($self->{_Coverage}{$id}{$pos}) }
	elsif (@_==4){ return($self->{_Coverage}{$id}{$pos}[$pop]) }
	elsif (@_==5){
		my $substring=substr($self->{_Coverage}{$id}{$pos}[$pop], $ind, 1);
		unless (defined $substring){ warn "Invalid access: $id, $pos, $pop, $ind!\n" }
		return( ord($substring) );
		}
	else { die "Invalid number of arguments!\n" }
	}


sub getGroups{
	my ($self, $group)=@_;
	if (@_==1){ return($self->{_Groups}) }
	elsif (@_==2){ return($self->{_Groups}[$group]) }
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

sub getIndnames{
	my ($self, $pop, $ind)=@_;
	if (@_==1){ return($self->{_Indnames}) }
	elsif (@_==3){ return($self->{_Indnames}[$pop][$ind]) }
	else { die "Wrong number of arguments!\n" }
	}


sub readvcfList{
	my ($self, $ref_loci, $mingtq, $suballeles, $recordbases, $recordcoverage, $onlyphased, $randomize, $verbose)=@_;

    # reset data:
    $self->{_Genotypes}=undef;
    $self->{_Bases}=undef;
    $self->{_Alleles}=undef;
    $self->{_Coverage}=undef;

	unless (defined $self->{_Files} && defined $self->{_Samples} ){ die "List and folder of vcf files have not been defined!\n" }
	my @popsizes=@{ $self->{_Popsizes} };
    my @filesizes=map { scalar(@$_) } (@{ $self->{_Samples} });

	my (@line, $pos, $iterator, @subinds_perpop, @reqphasing);
	if (defined $onlyphased && ref($onlyphased) eq 'ARRAY'){
		die "Wrong number of elements in onlyphased array!\n" unless (@$onlyphased==@{ $self->{_Files} });
		@reqphasing=@$onlyphased;
		}
	elsif (defined $onlyphased){ @reqphasing=map { $onlyphased } @{ $self->{_Files} } }
	else { @reqphasing=(0) x @{ $self->{_Files} } }

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if ( $ref_loci->getExcluded($chrom, $locus) );
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			if ($suballeles){
				$iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus);
				for my $ref_pop (@$iterator){
					push @subinds_perpop, { map {$_=>1} @$ref_pop };
					} 
				}
			FILE: for my $file (0..$#{ $self->{_Files} }){
				my $phased=$reqphasing[$file];
			    my @samples=@{ $self->{_Samples}[$file] };
				my @popmap=@{ $self->{_Popindices}[$file] };
				unless (@popmap==@samples){ die "Wrong length of population assignment map for vcf file $self->{_Files}[$file]!\n$!\n" }
				unless (@samples){
				    warn "No samples to process for file $file!\n";
				    next FILE;    # skip file if no samples to process.
				    }
				open my $vcffile, "-|", "tabix $self->{_Files}[$file] $chrom:$start-$end" or die "Could not open vcf file $self->{_Files}[$file]!\n$!\n";
				print STDERR "Aquiring genotypes for $chrom:$start-$end, vcf file $file.\n" if $verbose;
				LINE: while (<$vcffile>){
					if(/^[\w\-:]+\s/){
						@line=split(/\s/, $_);
						$pos=$line[1]-1;	# The sam and vcf formats are 1-based. All internal arrays and hashes use a 0-based indexing.
#						next LINE if ($self->{_Genotypes}{$line[0]}{$pos}[$file]);	# skip segregating position if already read in before.
                        if ($recordcoverage){
                            for my $pop (0..$#popsizes){
    						    unless (defined $self->{_Coverage}{$line[0]}{$pos}[$pop]){ $self->{_Coverage}{$line[0]}{$pos}[$pop]=(chr(0) x $popsizes[$pop]) }
    						    }
    						}
						if ($recordbases){
							my $basestring="$line[3]" . join('', split(',' , $line[4]) );
							my $current=$self->{_Bases}{$line[0]}{$pos};
							if (!defined $current){ $self->{_Bases}{$line[0]}{$pos}=$basestring }
							elsif (defined $current && $current ne $basestring){
								warn "$chrom:$pos, vcf file $file, bases $basestring don't match with bases from other populations $current!\n";
								my ($newbasestring, $ref_translation)=Misc::harmonizeBases($current, $basestring);
								$self->{_Bases}{$line[0]}{$pos}=$newbasestring;
								warn "Merged the two basestrings to $newbasestring.\n";
								warn "Applying the following translation of alleles to vcf file $file:\n";
								print STDERR "$_ to $ref_translation->{$_}\n" foreach (keys %$ref_translation);
#								print STDERR "Old line:\n", join("\t", @line), "\n";
								$line[4]=join(',', split('', substr($newbasestring, 1) ) );
								for my $call (@line[9..$#line]){
									$call=~s/([0-3])(\/|\|)([0-3])/$ref_translation->{$1}$2$ref_translation->{$3}/;
									}
#								print STDERR "New line:\n", join("\t", @line), "\n";
								}
							}
#							$self->{_Bases}{$line[0]}{$pos}="$line[3]" . join('', split(',' , $line[4]) ) }
						$self->{_Alleles}{$line[0]}{$pos}=0 unless defined $self->{_Alleles}{$line[0]}{$pos};
						my @tags=split(':', $line[8]);
						my %format=map {$tags[$_] => $_} (0..$#tags);
						my @index_perpop=(0) x @popsizes;	# record the current index of each population.
						for my $ind (0..$#samples){
						    my $pop=$popmap[$ind];
							my $col=$samples[$ind]+8;
							if ($col>$#line){ die "Invalid sample index specified for vcf file $self->{_Files}[$file]!\n" }
							if ($line[$col]=~/([0-3])(\/|\|)([0-3])/){
								my @call=split(':', $line[$col]);
								if ($phased && $2 ne '|' && $1 ne $3){
									if ($randomize){ $self->{_Genotypes}{$line[0]}{$pos}[$pop].=( int(rand(2)) ) ? $1 . $3 : $3 . $1 }
									else { $self->{_Genotypes}{$line[0]}{$pos}[$pop].=".." }
									}
								elsif ($mingtq && exists $format{'GQ'} && $call[$format{'GQ'}]<$mingtq){ $self->{_Genotypes}{$line[0]}{$pos}[$pop].=".." }
								else { $self->{_Genotypes}{$line[0]}{$pos}[$pop].=$1 . $3 }
								if ($recordcoverage && exists $format{'DP'}){
									my $coverage=$call[$format{'DP'}];
									substr($self->{_Coverage}{$line[0]}{$pos}[$pop], $index_perpop[$pop], 1)=chr($coverage<255 ? $coverage : 255);
									}
								if ($suballeles){	# only records allele states for subsampled individuals.
									if (exists $subinds_perpop[$pop]{ $index_perpop[$pop] }){ $self->setAlleles($line[0], $pos, $1, $3) }
									}
								else { $self->setAlleles($line[0], $pos, $1, $3) }
								}
							else { $self->{_Genotypes}{$line[0]}{$pos}[$pop].=".." }
							++$index_perpop[$pop];
							}
						}
					}
				close $vcffile;
				}
			}

		SITE: for my $seg ( keys %{ $self->{_Genotypes}{$chrom} } ){	# check if genotype vcf files have synchronized polymorphic positions.
			my $stringref=$self->getValue($chrom, $seg);
			for my $pop (0..$#popsizes){
				my $nalleles=2*$popsizes[$pop];
				if (defined $stringref->[$pop]){
					unless (length($stringref->[$pop])==$nalleles){
						warn "Genotype string for $chrom:$seg - population $pop has wrong length (", length($stringref->[$pop]), " vs. $nalleles)!",
							"Setting all genotypes in population to missing.\n";
						$stringref->[$pop]='.' x $nalleles;
						}
					}
				else {
#					warn "Genotype string for $chrom:$seg - population $pop is missing! Setting all genotypes to reference allele.\n";
					$stringref->[$pop]='0' x $nalleles;
					$self->setAlleles($chrom, $seg, '0');	# make sure that reference allele is set as present for the site.
					if ($recordcoverage){ $self->{_Coverage}{$chrom}{$seg}[$pop]=chr(0) x $nalleles }	# the use of the record coverage function in combination with non-synchronyzed vcf files is dangerous!
					}
				}
			}

		}
	return;
	}


1;

