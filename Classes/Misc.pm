package Misc;

use strict;
use warnings;



sub sum (@){
	my $rt=0;
	for my $i (@_){$rt+=$i}
	return($rt);
	}


sub min{
	my $ref=shift @_;
	my ($min_value, $min_index);
	my $index=0;
	for my $i (@$ref){
		if (defined $i){
			if (defined $min_value){ if ($i<$min_value){ $min_value=$i; $min_index=$index } }
			else { $min_value=$i; $min_index=$index }
			}
		++$index;
		}
	unless (defined $min_value){ warn "No valid elements in array!\n" }
	return($min_value, $min_index);
	}

sub minmax{
	my $ref=shift @_;
	my ($min, $max);
	for my $i (@$ref){
		if (defined $i){
			if (defined $min){ $min=$i if ($i<$min); $max=$i if ($i>$max) }
			else { $min=$i; $max=$i }
			}
		}
	unless (defined $min){ warn "No valid elements in array!\n" }
	return($min, $max);
	}


sub meansd{
	my $ref=shift @_;
	my ($n, $rt, $sum2)=(0, 0, 0);
	for my $i (@$ref){ if (defined $i){ $rt+=$i; ++$n } }
	unless ($n){ warn "No valid elements in array!\n"; return(undef, undef) }
	my $am=$rt / $n;
	for my $i (@$ref){ if (defined $i){ $sum2+=($i-$am)**2 } }
	my $sd=sqrt( $sum2 / $n );
	return($am, $sd);
	}


sub shuffle (@){
	my @a=\(@_);
	my $n;
	my $i=@_;
	map {
		$n=rand($i--);
		(${$a[$n]}, $a[$n]=$a[$i])[0];
		} @_;
	}


sub samplewithoutRep{
	my ($ref_array, $samplesize)=@_;
	unless (@_==2){ warn "Wrong number of arguments!\n"; return }
	my $i=@$ref_array;
	if ($samplesize>$i){ warn "Sample size ($samplesize) larger than number of elements in array ($i), setting samplesize to $i.\n"; $samplesize=$i }
	if ($samplesize==$i){ return(\@$ref_array) }
	my @a=\(@$ref_array);
	my $n;
	my @sample=map {
		$n=rand($i--);
		(${$a[$n]}, $a[$n]=$a[$i])[0];
		} (0..$samplesize-1);
	return(\@sample);
	}


sub expand{
	my $file=shift;
	$file=~s{(\d+)}{sprintf "%.8d", $1}eg; # expand all numbers to 8 digits
	return $file;
	}


sub countbits{
	unless (@_==2){die "Wrong number of arguments!\n"}
	my ($value, $nbits)=@_;
	my $count=0;
	my $lastbit;
	for my $bit (0..$nbits-1){ if ($value & (1 << $bit)){++$count; $lastbit=$bit} }
	return($count, $lastbit);
	}

sub countbits2{
	unless (@_==1){die "Wrong number of arguments!\n"}
	my ($value)=@_;
	my $count=0;
	my $bit=0;
	my $lastbit;
	while ($value){ if ($value & 1){ ++$count; $lastbit=$bit }; $value>>=1; ++$bit }
	return($count, $lastbit);
	}

sub getbits{
	unless (@_==1){die "Wrong number of arguments!\n"}
	my ($value)=@_;
	my @bits;
	my $bit=1;
	while ($value){ push @bits, $bit if ($value & 1); $value>>=1; ++$bit }
	return(@bits);
	}


sub openvcflist{
	my ($vcflist)=@_;
	unless (@_==1){ die "Wrong number of arguments!\n" }
	my (@vcffiles, @samples);
	open my $fh, "<", $vcflist or die "Could not open list of vcf files $vcflist!\n";
	while (<$fh>){
		if (/^#/){ next }
		chomp;
		my @line=split(/\s+/);
		push @vcffiles, $line[0];
		my @list;
		for my $entry ( split(',', $line[1]) ){
			my ($start, $end)=split('-', $entry);
			if (defined $end){ push @list, ($start..$end) }
			else { push @list, $start }
			}
		push @samples, \@list;
		}
	close $fh;
	return(\@vcffiles, \@samples);
	}


sub parseArglist{
	unless (@_==1){ die "Wrong number of arguments!\n" }
	my $argument=shift;
	my @list;
	for my $entry ( split(',', $argument) ){
		my ($start, $end)=split('-', $entry);
		if (defined $end){
			my ($sprefix, $sindex)=$start=~/([A-Za-z]*)(\d+)/;
			my ($eprefix, $eindex)=$end=~/([A-Za-z]*)(\d+)/;
			if ($eprefix && $eprefix ne $sprefix){ die "Unable to process argument list range $start-$end!" }
			for my $index ($sindex..$eindex){ push @list, "$sprefix" . "$index" }
			}
		else { push @list, $start }
		}
	print STDERR join(', ', @list), "\n";
	return(@list);
	}


sub readIndlist{
	unless (@_==1){ die "Wrong number of arguments!\n" }
	my $infile=shift;
	my @sampleids;
	open my $fh, "<", $infile or die "Could not open list of samples $infile!\n";
	while (<$fh>){
		next if (/^#/ || /^$/);
		chomp;
		my @line=split(/\s+/);
		push @sampleids, $line[0];
		}
	close $fh;
	return(@sampleids);
	}

sub readPoplist{
	my ($infile, $ref_poplabels)=@_;
	unless (@_>=1){ die "Wrong number of arguments!\n" }

	my (%poplabel_by_id, @sampleids, @poplabels, %popindices, @popsizes);
	if (defined $ref_poplabels){
	    @poplabels=@$ref_poplabels;
	    %popindices=map {$poplabels[$_]=>$_} (0..$#poplabels);
	    }

	open my $fh, "<", $infile or die "Could not open poplist $infile!\n";
	while (<$fh>){
		next if (/^#/ || /^$/);
		chomp;
		my @line=split(/\s+/);
		unless (@line>=2){ warn "Wrong number of entries for line $_\n!"; next }
		unless (defined $ref_poplabels || exists $popindices{$line[1]}){
		    push @poplabels, $line[1];
		    $popindices{$line[1]}=$#poplabels;
		    }
		push @sampleids, $line[0];
		unless (exists $popindices{$line[1]}){ die "Population $line[1] doesn't appear in list of population labels!\n" }
		$poplabel_by_id{$line[0]}=$line[1];
		++$popsizes[$popindices{$line[1]}];
		}
	close $fh;
	return(\%poplabel_by_id, \@sampleids, \@poplabels, \@popsizes);
	}

sub readGrouplist{
	unless (@_==1){ die "Wrong number of arguments!\n" }
	my $infile=shift;
	my (%grouplabel_by_poplabel, @grouplabels, %groupindices, %minpopind);

	open my $fh, "<", $infile or die "Could not open grouplist $infile!\n";
	while (<$fh>){
		next if (/^#/ || /^$/);
		chomp;
		my @line=split(/\s+/);
		unless (@line>=2){ warn "Wrong number of entries for line $_\n!"; next }
		unless (exists $groupindices{$line[1]}){
		    push @grouplabels, $line[1];
		    $groupindices{$line[1]}=$#grouplabels;
		    }
		$grouplabel_by_poplabel{$line[0]}=$line[1];
		$minpopind{$line[0]}=$line[2] if (defined $line[2]);
		}
	close $fh;
	return(\%grouplabel_by_poplabel, \@grouplabels, \%minpopind);
	}


sub getPopsums{
	my ($ref_popsizes)=@_;
	unless (@_==1){ die "Wrong number of arguments!\n" }
    my @popsums;
	my $rt=0;
	for my $popsize (@$ref_popsizes){
		push @popsums, $rt;
		$rt+=$popsize;
		}
	return(\@popsums);
	}


sub getPopindices{
	my ($ref_indnames, $ref_poplabels, $ref_popmap)=@_;
	if (@_==1){
    	my @fileindices;
    	for my $file (0..$#{ $ref_indnames }){
            push @fileindices, [ map { $file } @{ $ref_indnames->[$file] } ];
            }
        return(\@fileindices);	    
	    }
	elsif (@_==3){
        my %popidx=map { $ref_poplabels->[$_]=>$_ } (0..$#{ $ref_poplabels });
    	my @popindices;
    	for my $ref_file (@$ref_indnames){
            push @popindices, [ map { $popidx{ $ref_popmap->{$_} } } @$ref_file ];
            }
        return(\@popindices);
        }
    else { die "Wrong number of arguments!\n" }
    }

sub getPopfiles{
	unless (@_==1){ die "Wrong number of arguments!\n" }
	my $ref_popindices=shift;
	my %popfiles;
	for my $fileidx (0..$#{ $ref_popindices }){
		my %popindices=map { $_=>1 } @{ $ref_popindices->[$fileidx] };
		for my $popidx (keys %popindices){ 
			if (exists $popfiles{$popidx}){ warn "WARNING: Population index $popidx occurs in multiple files!\n" }
			$popfiles{$popidx}=$fileidx;
			}
		}
	return(\%popfiles);
	}

sub getAllindMap{
	my ($ref_popsizes, $ref_popmap_byfile, $verbose)=@_;
	unless (@_>=2){ die "Wrong number of arguments!\n" }
	
    my @allindmap;
    my $ref_popsums=getPopsums($ref_popsizes);

	for my $ref_popmap (@$ref_popmap_byfile){
        push @allindmap, [ map {$ref_popsums->[$_]++} @$ref_popmap ];
        if ($verbose){ print STDERR "Allindmap: ", join(',', @{ $allindmap[-1] }), "\n" }
        }
    return(\@allindmap);
    }


sub getIndnames{
	my ($vcflist, $coloffset, $ref_popmap)=@_;
	unless (@_>=1){ die "Wrong number of arguments!\n" }
	$coloffset//=8;	# in vcf files, the first sample name is at the 10th column.
	my ($ref_files, $ref_samples)=openvcflist($vcflist);
	my ($ref_indnames, $ref_popassignments, $ref_samples_updated)=readHeader($ref_files, $ref_samples, $coloffset, $ref_popmap);
	return($ref_indnames, $ref_files, $ref_samples_updated, $ref_popassignments);
	}

sub readHeader{
	my ($ref_files, $ref_samples, $coloffset, $ref_popmap)=@_;
	unless (@_>=2){ die "Wrong number of arguments!\n" }
	$coloffset//=8;	# in vcf files, the first sample name is at the 10th column.
	my (@indnames, @popassignments, @samples_by_file);
	for my $file (0..$#{$ref_files}){
		open my $vcfheader, "-|", "tabix $ref_files->[$file] -H" or die "Could not open vcf header $ref_files->[$file]!\n$!\n";
		VCF: while (<$vcfheader>){
			if(/^#CHROM/){
				chomp;
				my @header=split("\t", $_);
				unless ($#header>$coloffset){ die "VCF file $ref_files->[$file] doesn't contain any samples!\n" }
				my @samples=(!defined $ref_samples->[$file][0] || $ref_samples->[$file][0] eq 'all') ? (1..($#header-$coloffset)) : @{ $ref_samples->[$file] };
				for my $sample (@samples){
				    my $sampleid=$header[$coloffset+$sample];
				    my $popindex=(defined $ref_popmap) ? $ref_popmap->{$sampleid} : $file;
				    push @{ $indnames[$popindex] }, $sampleid;
				    push @{ $popassignments[$file] }, $popindex;
				    }
				push @samples_by_file, \@samples;
				last VCF;
				}
			}
		close $vcfheader;
		}
	return(\@indnames, \@popassignments, \@samples_by_file);
	}


sub getIndicesfromIndnames{
	my ($ref_allindnames, $ref_samples, $ref_subindnames, $verbose)=@_;
	unless (@_>=3){ die "Wrong number of arguments!\n" }

	my $nfiles=@$ref_allindnames;
    if ($ref_subindnames->[0] eq 'all'){
        return($ref_samples, $ref_allindnames);
        }
    	
	my %indices;
	for my $file (0..$nfiles-1){
		unless (@{ $ref_allindnames->[$file] }==@{ $ref_samples->[$file] } ){ die "Unequal sizes of indnames and samples arrays for file $file!\n" }
		$indices{ $ref_allindnames->[$file][$_] }=[ $file, $ref_samples->[$file][$_] ] foreach (0..$#{ $ref_allindnames->[$file] });
		}

	my @subsamples=map { [] } (1..$nfiles);
	my @subindnames=map { [] } (1..$nfiles);
	for my $indlabel (@$ref_subindnames){
		unless ( exists $indices{$indlabel} ){ warn "WARNING: Sublabel $indlabel does not exist in vcf files!\n"; next }
		my ($file, $index)=@{ $indices{$indlabel} };
		push @{ $subsamples[$file] }, $index;
		push @{ $subindnames[$file] }, $indlabel;		
		}
	
	if ($verbose){
    	print STDERR "$_: ", join( ',', @{ $subsamples[$_] } ), "\n" foreach (0..$#subsamples);
    	print STDERR "$_: ", join( ',', @{ $subindnames[$_] } ), "\n" foreach (0..$#subindnames);
        }
    
	return(\@subsamples, \@subindnames);
	}


sub getSubsamples{
	my ($ref_poplabels, $ref_indnames, $ref_samples, $ref_subpoplabels, $ref_selindnames, $ref_popmap, $ref_groupmap, $reorder)=@_;
	unless (@_>=6){ die "Wrong number of arguments!\n" }

    if ($ref_subpoplabels->[0] eq 'all' && $ref_selindnames->[0] eq 'all'){ return($ref_samples, $ref_indnames) }
    if ($ref_subpoplabels->[0] eq 'all'){ $ref_subpoplabels=$ref_poplabels }
    if ($ref_selindnames->[0] eq 'all'){
        my @allind;
        push @allind, @$_ foreach (@$ref_indnames);
        $ref_selindnames=\@allind;
        }
    
    my %poplabel_lookup=map { $_=>1 } (@$ref_subpoplabels);
    my %indname_lookup=map { $_=>1 } (@$ref_selindnames);

    my $nfiles=@$ref_indnames;
    my @subsamples=map { [] } (1..$nfiles);
    my @subindnames=map { [] } (1..$nfiles);

	if ($reorder){	# reorder individuals according to list of selected individuals:
		my %indices;
		for my $file (0..$nfiles-1){
			unless (@{ $ref_indnames->[$file] }==@{ $ref_samples->[$file] } ){ die "Unequal sizes of indnames and samples arrays for file $file!\n" }
			$indices{ $ref_indnames->[$file][$_] }=[ $file, $ref_samples->[$file][$_] ] foreach (0..$#{ $ref_indnames->[$file] });
			}

		for my $indname (@$ref_selindnames){
			unless (exists $indices{$indname} ){ warn "Sublabel $indname does not exist in poplist!\n"; next }
            my $poplabel=$ref_popmap->{$indname};
            my $grouplabel=(defined $ref_groupmap && exists $ref_groupmap->{$poplabel}) ? $ref_groupmap->{$poplabel} : $poplabel;
			if ($poplabel_lookup{$poplabel} || $poplabel_lookup{$grouplabel}){
				my ($file, $sample)=@{ $indices{$indname} };
				push @{ $subsamples[$file] }, $sample;
				push @{ $subindnames[$file] }, $indname;
				}
			}
		}
	else {	# loop through all individuals to conserve order of samples according to poplist:
	    for my $file (0..$nfiles-1){
 	    	for my $ind (0..$#{ $ref_samples->[$file] }){
            	my $indname=$ref_indnames->[$file][$ind];
            	my $poplabel=$ref_popmap->{$indname};
            	my $grouplabel=(defined $ref_groupmap && exists $ref_groupmap->{$poplabel}) ? $ref_groupmap->{$poplabel} : $poplabel;
    	    	if (($poplabel_lookup{$poplabel} || $poplabel_lookup{$grouplabel}) && $indname_lookup{$indname}){
                	push @{ $subsamples[$file] }, $ref_samples->[$file][$ind];
                	push @{ $subindnames[$file] }, $indname;
                	}
    	        }
    	    }   
        } 
    return(\@subsamples, \@subindnames);
    }


sub getPops{
	my ($ref_indnames, $ref_poplabels, $ref_popmap, $ref_grouplabels, $ref_groupmap, $verbose)=@_;
	unless (@_>=3){ die "Wrong number of arguments!\n" }	
	
    unless (defined $ref_grouplabels && defined $ref_groupmap){
        $ref_grouplabels=$ref_poplabels;
        $ref_groupmap={ map {$_=>$_} @$ref_poplabels };
        }
    
	# count number of individuals per population:
    my %indnames_bypop;
	for my $ref_file (@$ref_indnames){
	    for my $indname (@$ref_file){
    	    push @{ $indnames_bypop{ $ref_popmap->{$indname} } }, $indname;
	        }
	    }
	
	# modify population label array:
    my (@subpoplabels, @subpopsizes, @subindnames, %poplabels_bygroup);
	for my $poplabel (@$ref_poplabels){
	    if (exists $indnames_bypop{$poplabel}){
			push @subpoplabels, $poplabel;
			push @subpopsizes, scalar(@{ $indnames_bypop{$poplabel} });
			push @subindnames, $indnames_bypop{$poplabel};
			push @{ $poplabels_bygroup{ $ref_groupmap->{$poplabel} } }, $poplabel if (exists $ref_groupmap->{$poplabel});
			}
	    }

    # build population assignment array:
    my $ref_subpopassignments=getPopindices($ref_indnames, \@subpoplabels, $ref_popmap);
    my %subpopidx=map { $subpoplabels[$_]=>$_ } (0..$#subpoplabels);

	# modify group label array:
    my @subgrouplabels=grep {exists $poplabels_bygroup{$_} } @$ref_grouplabels;

    # build group assignment array:
	my @groups;
	for my $grouplabel (@subgrouplabels){
        push @groups, [ map { $subpopidx{$_} } @{ $poplabels_bygroup{$grouplabel} } ];
        print STDERR "Group $grouplabel: ", join(',', @{ $groups[-1] }), "\n" if $verbose;
        }
    
    if ($verbose){
        print STDERR "File\tIndividualID\tGrouplabel\tPoplabel\tOldPopidx\tNewPopidx\n";
        my %popidx=map { $ref_poplabels->[$_]=>$_ } (0..$#{ $ref_poplabels });
        for my $file (0..$#{ $ref_indnames }){
            for my $ind (0..$#{ $ref_indnames->[$file] }){
                my $indname=$ref_indnames->[$file][$ind];
                my $poplabel=$ref_popmap->{$indname};
                my $grouplabel=(defined $ref_groupmap && exists $ref_groupmap->{$poplabel}) ? $ref_groupmap->{$poplabel} : $poplabel;
                print STDERR "$file\t$indname\t$grouplabel\t$poplabel\t$popidx{$poplabel}\t$ref_subpopassignments->[$file][$ind]\n";
                }
            }
        }

    return(\@subpopsizes, \@subpoplabels, \@groups, \@subgrouplabels, \@subindnames, $ref_subpopassignments);
    }


sub iupac{
	my ($base1, $base2)=@_;
	unless (@_==2){die "Wrong number of arguments!\n"}

	if ($base1 eq "A"){
		if ($base2 eq "T"){return "W"}
		if ($base2 eq "G"){return "R"}
		if ($base2 eq "C"){return "M"}
		if ($base2 eq "A"){return "A"}
		}
	if ($base1 eq "T"){
		if ($base2 eq "A"){return "W"}
		if ($base2 eq "G"){return "K"}
		if ($base2 eq "C"){return "Y"}
		if ($base2 eq "T"){return "T"}
		}
	if ($base1 eq "G"){
		if ($base2 eq "A"){return "R"}
		if ($base2 eq "T"){return "K"}
		if ($base2 eq "C"){return "S"}
		if ($base2 eq "G"){return "G"}
		}
	if ($base1 eq "C"){
		if ($base2 eq "A"){return "M"}
		if ($base2 eq "T"){return "Y"}
		if ($base2 eq "G"){return "S"}
		if ($base2 eq "C"){return "C"}
		}
	warn "Input bases not A,T,G or C! $base1, $base2\n";
	return("N");
	}


sub getIupacConsensus{
	my ($ref_seq1, $ref_seq2)=@_;
	unless (@_==2){die "Wrong number of arguments!\n"}
	
	# make hardcopies of sequences:
	my $seq1=$$ref_seq1;
	my $seq2=$$ref_seq2;
	unless (length($seq1)==length($seq2)){ die "Unequal length of unphased haplotype strings!\n" }

	my $consensus;
	while (length($seq1)>0){
		my ($base1, $base2)=(substr($seq1, 0, 1, ""), substr($seq2, 0, 1, "") );	# extract first bases in strings and chew them up for faster access.
		if ($base1 eq $base2){	# if both bases identical, take the first one for both haplotype strings.
			$consensus.=$base1;
			}
		else {
			$consensus.=iupac($base1, $base2);	# randomly select one of the two bases for the first haplotype string.
			}
		}
	unless (length($consensus)==length($$ref_seq1) && length($consensus)==length($$ref_seq2)){ die "Unequal length of input and output sequence strings!\n" }
	return(\$consensus);
	}


sub getCodonRef{
	my(%codon_hash)=(
		'TCA' => 'S', # Serine
		'TCC' => 'S', # Serine
		'TCG' => 'S', # Serine
		'TCT' => 'S', # Serine
		'TTC' => 'F', # Phenylalanine
		'TTT' => 'F', # Phenylalanine
		'TTA' => 'L', # Leucine
		'TTG' => 'L', # Leucine
		'TAC' => 'Y', # Tyrosine
		'TAT' => 'Y', # Tyrosine
		'TAA' => '_', # Stop
		'TAG' => '_', # Stop
		'TGC' => 'C', # Cysteine
		'TGT' => 'C', # Cysteine
		'TGA' => '_', # Stop
		'TGG' => 'W', # Tryptophan
		'CTA' => 'L', # Leucine
		'CTC' => 'L', # Leucine
		'CTG' => 'L', # Leucine
		'CTT' => 'L', # Leucine
		'CCA' => 'P', # Proline
		'CAT' => 'H', # Histidine
		'CAA' => 'Q', # Glutamine
		'CAG' => 'Q', # Glutamine
		'CGA' => 'R', # Arginine
		'CGC' => 'R', # Arginine
		'CGG' => 'R', # Arginine
		'CGT' => 'R', # Arginine
		'ATA' => 'I', # Isoleucine
		'ATC' => 'I', # Isoleucine
		'ATT' => 'I', # Isoleucine
		'ATG' => 'M', # Methionine
		'ACA' => 'T', # Threonine
		'ACC' => 'T', # Threonine
		'ACG' => 'T', # Threonine
		'ACT' => 'T', # Threonine
		'AAC' => 'N', # Asparagine
		'AAT' => 'N', # Asparagine
		'AAA' => 'K', # Lysine
		'AAG' => 'K', # Lysine
		'AGC' => 'S', # Serine
		'AGT' => 'S', # Serine
		'AGA' => 'R', # Arginine
		'AGG' => 'R', # Arginine
		'CCC' => 'P', # Proline
		'CCG' => 'P', # Proline
		'CCT' => 'P', # Proline
		'CAC' => 'H', # Histidine
		'GTA' => 'V', # Valine
		'GTC' => 'V', # Valine
		'GTG' => 'V', # Valine
		'GTT' => 'V', # Valine
		'GCA' => 'A', # Alanine
		'GCC' => 'A', # Alanine
		'GCG' => 'A', # Alanine
		'GCT' => 'A', # Alanine
		'GAC' => 'D', # Aspartic Acid
		'GAT' => 'D', # Aspartic Acid
		'GAA' => 'E', # Glutamic Acid
		'GAG' => 'E', # Glutamic Acid
		'GGA' => 'G', # Glycine
		'GGC' => 'G', # Glycine
		'GGG' => 'G', # Glycine
		'GGT' => 'G'  # Glycine
		);
	return(\%codon_hash);
	}


sub getDegRef{
	my $ref_codon_hash=getCodonRef();
	my %degenerancy;

	while( my ($codon, $aa)=each %$ref_codon_hash ){
		if ($aa eq '_'){ $degenerancy{$codon}='stop'; next } 
		my @deg=(0,0,0);
		for my $pos (0..2){
			for my $base ('A','T','C','G'){
				my $altcodon=$codon;
				substr($altcodon, $pos, 1)=$base;
				my $altaa=$ref_codon_hash->{$altcodon};
#				print STDERR "$altcodon - $altaa\n";
				++$deg[$pos] if ($aa eq $altaa);
				}
			die "Error in generating degeneracy hash!\n" unless ($deg[$pos]);
			}
		$degenerancy{$codon}=join('', @deg);
#		print STDERR "codon: $codon - degenerancy: ", $degenerancy{$codon}, "\n";
		}

	return(\%degenerancy);
	}


sub reverseComplement{	# IUPAC characters will not be properly reverse complemented!
	my ($ref_sequence)=@_;

	my $revseq=reverse($$ref_sequence);
	my $nrep=($revseq=~tr/AaTtGgCcN/TtAaCcGgN/);
	$nrep+=($revseq=~tr/WwRrMmKkYySs/N/);
	unless ($nrep==length($$ref_sequence) ){ warn "Invalid characters found in sequence!\n" }
	unless (length($revseq)==length($$ref_sequence) ){ die "Length of sequence after reverse complementing not the same as length of input sequence!\n" }

	return(\$revseq);
	}


sub getWeakStrong{
	my ($ancbase, $derbase)=@_;
	unless (@_==2){die "Wrong number of arguments!\n"}
	my $cat;

	if ($ancbase eq 'A' || $ancbase eq 'T'){
		if ($derbase eq 'A' || $derbase eq 'T'){ $cat=0 }
		elsif ($derbase eq 'G' || $derbase eq 'C'){ $cat=1 }
		}
	elsif ($ancbase eq 'G' || $ancbase eq 'C'){
		if ($derbase eq 'A' || $derbase eq 'T'){ $cat=2 }
		elsif ($derbase eq 'G' || $derbase eq 'C'){ $cat=3 }
		}

	return($cat);
	}


sub getCpGString{
	my ($ref_seqstring, $genotypes, $chrom, $start, $end)=@_;

	my ($cpg, $nextc);
	my $seqstring=$$ref_seqstring;	# copy string from reference.
	my $base=uc( substr($seqstring, 0, 1, "") );	# chew up first base.
	my $prevc=1;	# Discard first base because n-1 state is missing.
	POSITION: for my $pos ($start..$end-1){
		my $offset=$pos-$start;
		my $nowc=0; my $nextg=0;
		if (defined $nextc){ $nowc=$nextc }
		else {
			if ($base eq 'C'){ $nowc=1 }
			elsif ( defined $genotypes->getValue($chrom, $pos) ){
				if ( $genotypes->NAlleles($chrom, $pos)>2 ){ $nowc=1 }
				elsif ( uc( $genotypes->getBase($chrom, $pos, 1) ) eq 'C'){ $nowc=1 };
				}
			}
		$nextc=0;
		my $nextbase=uc( substr($seqstring, 0, 1, "") );	# chew up next base of string.
		if ($nextbase eq 'G'){ $nextg=1 }
		elsif ($nextbase eq 'C'){ $nextc=1 }
		elsif ( defined $genotypes->getValue($chrom, $pos+1) ){
			if ( $genotypes->NAlleles($chrom, $pos+1)>2 ){ $nextg=1; $nextc=1; warn "More than two alleles at position $chrom:$pos!\n" }	# non-biallelic bases are treated like containing C and G.
			else {
				my $nextalt=uc( $genotypes->getBase($chrom, $pos+1, 1) );
				if ($nextalt eq 'G'){ $nextg=1 }
				elsif ($nextalt eq 'C'){ $nextc=1 }
				}
			}
		if ($prevc || $nextg){ $cpg.='0' }
		else { $cpg.='1' }
		$prevc=$nowc;
		$base=$nextbase;
		}
	$cpg.='0';	# Discard last base because n+1 state is missing.
#	print STDERR "$chrom:$start-$end\n$cpg\n";
	return(\$cpg);
	}


sub maskInterval{
	my ($start, $end, @ref_maskarray)=@_;

	my @combinedmask;
	push @combinedmask, @{$_} foreach @ref_maskarray;
	my @positions;
	my $tempstart=$start;
	my $complete=0;
	REG: for my $reg ( sort { $a->{'start'}<=>$b->{'start'} } @combinedmask ){
		if ( $reg->{'end'}<$tempstart ){ next REG }
		elsif ( $reg->{'start'}>$end ){ last REG }
		elsif ( $reg->{'start'}>$tempstart ){
			push @positions, ( $tempstart..$reg->{'start'}-1 );
			$tempstart=$reg->{'end'}+1;
			if ($tempstart>$end){ $complete=1; last REG }
			}
		elsif ( $reg->{'start'}<=$tempstart ){
			$tempstart=$reg->{'end'}+1 if $reg->{'end'}>=$tempstart;
			if ($tempstart>$end){ $complete=1; last REG }
			}
		}
	push @positions, ($tempstart..$end) unless ($complete);	# return full range if all mask regions are smaller than end of range.
	return(\@positions);
	}

sub maskIntervalRanges{
	my ($start, $end, @ref_maskarray)=@_;

	my @combinedmask;
	push @combinedmask, @{$_} foreach @ref_maskarray;
	my @ranges;
	my $tempstart=$start;
	my $complete=0;
	REG: for my $reg ( sort { $a->{'start'}<=>$b->{'start'} } @combinedmask ){
		if ( $reg->{'end'}<$tempstart ){ next REG }
		elsif ( $reg->{'start'}>$end ){ last REG }
		elsif ( $reg->{'start'}>$tempstart ){
			push @ranges, [$tempstart, $reg->{'start'}-1];
			$tempstart=$reg->{'end'}+1;
			if ($tempstart>$end){ $complete=1; last REG }
			}
		elsif ( $reg->{'start'}<=$tempstart ){
			$tempstart=$reg->{'end'}+1 if $reg->{'end'}>=$tempstart;
			if ($tempstart>$end){ $complete=1; last REG }
			}
		}
	push @ranges, [$tempstart, $end] unless ($complete);	# return full range if all mask regions are smaller than end of range.
	return(\@ranges);
	}


sub printSFSforDFEalpha{
	my ($genes, $introns, $outfile, $folded)=@_;
	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", $outfile or die "Could not open output file!\n"}

	my @nalleles;
	print $output scalar @{ $genes->{_SFS} }, "\n", scalar @{ $genes->{_SFS}[0] }, "\n";
	for my $set (0..$#{ $genes->{_SFS}[0] } ){
		push @nalleles, $folded ? (scalar @{ $genes->{_SFS}[0][$set] }-1)*2 : scalar @{ $genes->{_SFS}[0][$set] }-1;
		print $output "$nalleles[-1]\n";
		}

	for my $pop (0..$#{ $genes->{_SFS} } ){
		print $output $pop+1, "\n???\n???\n???\n???\n";
		for my $set (0..$#{ $genes->{_SFS}[$pop] } ){
			for my $nder (0..$nalleles[$set]){
				my $count=defined $genes->{_SFS}[$pop][$set][$nder] ? $genes->{_SFS}[$pop][$set][$nder] : 0;	# this needs to be the selected SFS.
				print $output "$count ";
				}
			print $output "\n";
			for my $nder (0..$nalleles[$set]){
				my $count=defined $introns->{_SFS}[$pop][$set][$nder] ? $introns->{_SFS}[$pop][$set][$nder] : 0;	# this needs to be the neutral SFS.
				print $output "$count ";
				}
			print $output "\n";
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}


sub setSubArrays{
	my ($ref_sublabels, $ref_poplabels, $ref_nindividuals, $ref_minind, $ref_popgroups, $ref_grouplabels, $ref_groupminind)=@_;
	my %sublabels=map { $_=>1 } (@$ref_sublabels);
	my (@subpopindices, @subnindividuals, @subminind, @subpopgroups, @subpoplabels, @subgrouplabels, @subgroupminind);
	my $newindex=0;
	for my $group ( 0..$#{ $ref_popgroups } ){
		my @subgroup;
		for my $pop (@{ $ref_popgroups->[$group] }){
			if ( $sublabels{ $ref_poplabels->[$pop] } || $sublabels{ $ref_grouplabels->[$group] } ){
				push @subgroup, $newindex++;
				push @subpoplabels, $ref_poplabels->[$pop];
				push @subpopindices, $pop;
				push @subnindividuals, $ref_nindividuals->[$pop];
				if ( ref($ref_minind->[0]) eq "ARRAY" ){ push @{ $subminind[$_] }, $ref_minind->[$_][$pop] foreach (0..$#{ $ref_minind }) }	# to support multiple sets of minimum individuals.
				else { push @subminind, $ref_minind->[$pop] }
				}
			}
		if (@subgroup){
			push @subpopgroups, \@subgroup;
			push @subgrouplabels, $ref_grouplabels->[$group];
			if (defined $ref_groupminind){
				push @subgroupminind, $ref_groupminind->[$group];
				if ($ref_groupminind->[$group]>sum(@subnindividuals[@subgroup]) ){ 
					die "Minimum number of individuals required in group $group is larger than total number of individuals in subset of populations! ($ref_groupminind->[$group] vs. ", sum(@subnindividuals[@subgroup]), ")\n";
					}
				}
			}
		}

	# overwrite original population arrays:
	@$ref_nindividuals=@subnindividuals;
	@$ref_minind=@subminind;
	@$ref_popgroups=@subpopgroups;
	@$ref_poplabels=@subpoplabels;
	@$ref_grouplabels=@subgrouplabels;
	@$ref_groupminind=@subgroupminind if (defined $ref_groupminind);

	return(@subpopindices);
	}


sub contract_linear_paths{
	my $tree=shift;
	my $reroot=shift;
	my @remove_internal_ids;
	foreach my $node ($tree->get_nodes){
		if ($node->ancestor && $node->each_Descendent==1){
			push(@remove_internal_ids, $node->internal_id);
			}
		}
	$tree->splice(-remove_internal_id => \@remove_internal_ids, -preserve_lengths => 1) if @remove_internal_ids;
	if ($reroot){
		my $root=$tree->get_root_node;
		my @descs=$root->each_Descendent;
		if (@descs==1){
			my $new_root=shift(@descs);
			$tree->set_root_node($new_root);
			$new_root->ancestor(undef);
			}
		}
	}


sub setMissing{
	unless (@_==1){ die "Wrong number of arguments!\n" }
	my ($ref_line)=@_;
	$ref_line->[7]="MLEAC=0;MLEAF=0.000";
	$ref_line->[$_]="./." foreach ( 9..$#{$ref_line} );
	return;
	}

sub setFilter{
	unless (@_==3){ die "Wrong number of arguments!\n" }
	my ($ref_line, $status, $ref_bitflags)=@_;
	my @bits=getbits($status);
	my $filterstring=join(';', map { $ref_bitflags->{$_} } @bits);
	$filterstring="PASS" unless ( length($filterstring)>0);
	if ($ref_line->[6] eq "PASS" || $ref_line->[6] eq "."){ $ref_line->[6]=$filterstring }
	else { $ref_line->[6].=";$filterstring" }
	return;
	}


sub subsampleVCF{
	unless (@_==3){ die "Wrong number of arguments!\n" }
	my ($ref_line, $ref_subind, $haplotize)=@_;
	my @line=@{ $ref_line }[0..8];
	my @genotypes;
	my @ac=(0) x 4;
	my $tc=0;
	for my $ind (map { @{ $ref_line }[$_+8] } @$ref_subind){
		if ($ind=~/([0-3])\/([0-3])/){
			++$ac[$1];
			++$ac[$2];
			$tc+=2;
			}
		push @genotypes, $ind;
		}
	my @af=map { sprintf("%.3f", $_/$tc) } @ac;
	my @info;

	if ($haplotize){
		my $major;
		my @majors;
		my $maxvalue=0;
		for my $index (0..3){
			if ($ac[$index]>$maxvalue){ @majors=($index) }
			elsif ($ac[$index]==$maxvalue){ push @majors, $index }
			}
		if (@majors>1){ $major=$majors[ int( rand(scalar(@majors)) ) ] }
		else { $major=$majors[0] }	
		push @line, "${major}/${major}";
		$info[0]=$major==0 ? "AC=0" : "AC=1";
		$info[1]=$major==0 ? "AF=0.000" : "AF=1.000";
		$info[2]="AN=1";
		}
	else {
		push @line, @genotypes;
		$info[0]="AC=" . join(',', @ac[1..3]);
		$info[1]="AF=" . join(',', @af[1..3]);
		$info[2]="AN=" . $tc;
		}
	$line[7]=join(';', @info);
	return(\@line);
	}


sub harmonizeBases{
	unless (@_==2){ die "Wrong number of arguments!\n" }
	my ($bases1, $bases2)=@_;

	my @bases1=split('', $bases1);
	my @bases2=split('', $bases2);
	die "First bases of the two basestrings don't match! Different reference genomes were used to generate vcf files!\n" unless ($bases1[0] eq $bases2[0]);

	my $newbasestring=$bases1;
	my %translation;
	my %map=map {$bases1[$_]=>$_} (0..$#bases1);

	my $nextidx=@bases1;
	for my $i (0..$#bases2){
		my $testbase=$bases2[$i];
		if (exists $map{$testbase}){
			if ($map{$testbase}==$i){ $translation{$i}=$i }
			else { $translation{$i}=$map{$testbase} }
			}
		else { $map{$testbase}=$nextidx; $newbasestring.=$testbase; $translation{$i}=$nextidx; ++$nextidx }
		}

	return($newbasestring, \%translation);
	}


1;

