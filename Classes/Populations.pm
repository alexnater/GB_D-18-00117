package Populations;

use strict;
use warnings;
use Classes::Misc;


sub new{
	my ($class, $poplist_file, $grouplist_file, $ref_subpops, $ref_subinds, $reorder)=@_;
	unless (@_>=2){ die "ERROR: Wrong number of arguments!\n" }
	my $self={};
	bless $self, $class;
	$self->readPoplist($poplist_file);

	if (defined $grouplist_file && $grouplist_file ne 'none'){
		$self->readGrouplist($grouplist_file);
		}
	if (defined $ref_subinds){
    	my @subindividuals=@$ref_subinds;
	    if ($ref_subinds->[0]=~/file=(.+)$/){
		    # read file with list of selected individuals:
            print STDERR "Reading list of individuals from file $1.\n";
            @subindividuals=Misc::readIndlist($1);
            }
        $self->setSubset($ref_subpops, \@subindividuals, $reorder);
        }
    elsif (defined $ref_subpops){
	    $self->setSubset($ref_subpops, ['all'], $reorder);
		}
	$self->setGroups();
	return $self;
	}


sub getSampleids{
	my $self=shift;
	return($self->{_Sampleids});
	}

sub getSampleids_bypop{
	my ($self, $pop)=@_;
	if (@_==1){ return($self->{_Sampleids_bypop}) }
	elsif (@_==2){ return($self->{_Sampleids_bypop}[$pop]) }
	else { die "ERROR: Wrong number of arguments!\n" }
	}

sub getPopsizes{
	my ($self, $popidx)=@_;
	if (@_==1){ return($self->{_Popsizes}) }
	elsif (@_==2){ return($self->{_Popsizes}[$popidx]) }
	else { die "ERROR: Wrong number of arguments!\n" }
	}

sub getPoplabels{
	my ($self, $popidx)=@_;
	if (@_==1){ return($self->{_Poplabels}) }
	elsif (@_==2){ return($self->{_Poplabels}[$popidx]) }
	else { die "ERROR: Wrong number of arguments!\n" }
	}

sub getPop{
	my ($self, $sampleid)=@_;
	unless (@_==2){ die "ERROR: Wrong number of arguments!\n" }
	my $poplabel=$self->{_Popmap}{$sampleid};
	return($poplabel, $self->{_Popidx}{$poplabel});
	}

sub getPopmap{
	my $self=shift;
	return($self->{_Popmap});
	}

sub getGrouplabels{
	my ($self, $groupidx)=@_;
	if (@_==1){ return($self->{_Grouplabels}) }
	elsif (@_==2){ return($self->{_Grouplabels}[$groupidx]) }
	else { die "ERROR: Wrong number of arguments!\n" }
	}

sub getGroup{
	my ($self, $poplabel)=@_;
	unless (@_==2){ die "ERROR: Wrong number of arguments!\n" }
	my $grouplabel=$self->{_Groupmap}{$poplabel};
	return($grouplabel, $self->{_Groupidx}{$grouplabel});
	}

sub getGroupmap{
	my $self=shift;
	return($self->{_Groupmap});
	}

sub getGroups{
	my $self=shift;
	return($self->{_Groups});
	}

sub getGroupsizes{
	my ($self, $groupidx)=@_;
	if (@_==1){ return($self->{_Groupsizes}) }
	elsif (@_==2){ return($self->{_Groupsizes}[$groupidx]) }
	else { die "ERROR: Wrong number of arguments!\n" }
	}

sub getSampleinfo{
	my ($self, $sampleid)=@_;
	if (@_==1){ return($self->{_Sampleinfo}) }
	elsif (@_==2){ return($self->{_Sampleinfo}{$sampleid}) }
	else { die "ERROR: Wrong number of arguments!\n" }
	}


sub readPoplist{
	my ($self, $infile, $ref_poplabels)=@_;
	unless (@_>=2){ die "ERROR: Wrong number of arguments!\n" }

	my (%poplabel_by_id, @sampleids, @sampleids_bypop, @poplabels, %popindices, @popsizes, %infomap);
	if (defined $ref_poplabels){
	    @poplabels=@$ref_poplabels;
	    %popindices=map {$poplabels[$_]=>$_} (0..$#poplabels);
	    }

	open my $fh, "<", $infile or die "ERROR: Could not open poplist $infile!\n";
	while (<$fh>){
		next if (/^#/ || /^$/);
		chomp;
		my @line=split(/\s+/);
		unless (@line>=2){ warn "Wrong number of entries for line $_\n!"; next }
		unless (defined $ref_poplabels || exists $popindices{$line[1]}){
		    push @poplabels, $line[1];
		    $popindices{$line[1]}=$#poplabels;
		    }
		unless (exists $popindices{$line[1]}){ die "ERROR: Population $line[1] doesn't appear in list of population labels!\n" }
		push @sampleids, $line[0];
		push @{ $sampleids_bypop[$popindices{$line[1]}] }, $line[0];
		$poplabel_by_id{$line[0]}=$line[1];
		++$popsizes[$popindices{$line[1]}];
		$infomap{$line[0]}=$line[2] if (@line>2);
		}
	close $fh;
	
	$self->changePoplabels(\@poplabels, \@popsizes, \@sampleids_bypop);
	$self->{_Sampleids}=\@sampleids;
	$self->{_Popmap}=\%poplabel_by_id;
	$self->{_Sampleinfo}=\%infomap;
	
	# set group data to populations:
	$self->changeGrouplabels(\@poplabels);
	$self->{_Groupmap}={ map {$_ => $_} @poplabels };
	return;
	}


sub readGrouplist{
	my ($self, $infile)=@_;
	unless (@_==2){ die "ERROR: Wrong number of arguments!\n" }

	my (%grouplabel_by_poplabel, @grouplabels, %groupindices, %minpopind);

	open my $fh, "<", $infile or die "ERROR: Could not open grouplist $infile!\n";
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

	# if population has no group assignment, set grouplabel to poplabel:
	for my $poplabel (@{ $self->{_Poplabels} }){
		unless (exists $grouplabel_by_poplabel{$poplabel}){
			$grouplabel_by_poplabel{$poplabel}=$poplabel;
			push @grouplabels, $poplabel;
			}
		}
	
	$self->changeGrouplabels(\@grouplabels);
	$self->{_Minpopsizes}=\%minpopind;
	$self->{_Groupmap}=\%grouplabel_by_poplabel;
	return;
	}


sub readIndlist{
	my ($self, $infile)=@_;
	unless (@_==2){ die "ERROR: Wrong number of arguments!\n" }

	my @sampleids;
	open my $fh, "<", $infile or die "ERROR: Could not open list of samples $infile!\n";
	while (<$fh>){
		next if (/^#/ || /^$/);
		chomp;
		my @line=split(/\s+/);
		push @sampleids, $line[0];
		}
	close $fh;
	
	$self->{_Subsampleids}=\@sampleids;
	return;
	}


sub setSubset{
	my ($self, $ref_subpoplabels, $ref_subindnames, $reorder)=@_;
	unless (@_>=3){ die "ERROR: Wrong number of arguments!\n" }

    return if ($ref_subpoplabels->[0] eq 'all' && $ref_subindnames->[0] eq 'all');
    if ($ref_subpoplabels->[0] eq 'all'){ $ref_subpoplabels=$self->{_Poplabels} }
    if ($ref_subindnames->[0] eq 'all'){ $ref_subindnames=$self->{_Sampleids} }
    
    my %poplabel_lookup=map { $_=>1 } (@$ref_subpoplabels);
    my (@subsampleids, %indnames_bypop);

	if ($reorder){	# reorder individuals within populations according to list of selected individuals:
		for my $indname (@$ref_subindnames){
			unless (exists $self->{_Popmap}{$indname}){ warn "WARNING: Sublabel $indname does not exist in poplist!\n"; next }
            my $poplabel=$self->{_Popmap}{$indname};
            my $grouplabel=$self->{_Groupmap}{$poplabel};
			if ($poplabel_lookup{$poplabel} || $poplabel_lookup{$grouplabel}){
		    	push @subsampleids, $indname;
				push @{ $indnames_bypop{$poplabel} }, $indname;
				}
			}
		}
	else {	# loop through all individuals to conserve order of samples according to poplist:
	    my %subindname_lookup=map { $_=>1 } (@$ref_subindnames);
 	    for my $indname (@{ $self->{_Sampleids} }){
           	my $poplabel=$self->{_Popmap}{$indname};
           	my $grouplabel=$self->{_Groupmap}{$poplabel};
   	    	if (($poplabel_lookup{$poplabel} || $poplabel_lookup{$grouplabel}) && $subindname_lookup{$indname}){
				push @subsampleids, $indname;
				push @{ $indnames_bypop{$poplabel} }, $indname; 
    	        }
    	    }   
        }
	
	# modify population label array:
    my (@subpoplabels, @subpopsizes, @subindnames_bypop);
	for my $poplabel (@{ $self->{_Poplabels} }){
	    if (exists $indnames_bypop{$poplabel}){
			push @subpoplabels, $poplabel;
			push @subpopsizes, scalar(@{ $indnames_bypop{$poplabel} });
			push @subindnames_bypop, $indnames_bypop{$poplabel};
			}
	    }

  	$self->changePoplabels(\@subpoplabels, \@subpopsizes, \@subindnames_bypop);
	$self->{_Sampleids}=\@subsampleids;
    return;
    }


sub SetSubset_byCoverage{
	my ($self, $ind_perpop)=@_;
	unless (@_==2){ die "ERROR: Wrong number of arguments!\n" }
	if ($ind_perpop eq "all"){ return }
	unless (defined $self->{_Sampleinfo}){ die "ERROR: No coverage information provided in poplist!\n" }
	my (@subindnames_bypop, @subsampleids, @popsizes);
	my $popidx=0;
	for my $ref_sampleids (@{ $self->{_Sampleids_bypop} }){
		if (@$ref_sampleids<$ind_perpop){
			warn "WARNING: Not enough individuals for population $popidx (", scalar(@$ref_sampleids), " vs. $ind_perpop)!\n";
			$ind_perpop=scalar(@$ref_sampleids);
			}
		my @popsubids=(sort { $self->{_Sampleinfo}{$b}<=>$self->{_Sampleinfo}{$a} } @$ref_sampleids)[0..$ind_perpop-1];
		print STDERR "Population $self->{_Poplabels}[$popidx]:\n";
		print STDERR "$_\t$self->{_Sampleinfo}{$_}\n" foreach (@popsubids);
		push @subindnames_bypop, \@popsubids;
		push @subsampleids, @popsubids;
		push @popsizes, scalar(@popsubids);
		++$popidx;
		}
	
	$self->{_Popsizes}=\@popsizes;
	$self->{_Sampleids}=\@subsampleids;
	$self->{_Sampleids_bypop}=\@subindnames_bypop;
	return;
	}


sub setGroups{
	my $self=shift;

	my (%poplabels_bygroup, %popindices_bygroup, @subgrouplabels, @popgroups);

	for my $poplabel (@{ $self->{_Poplabels} }){
		my $grouplabel=$self->{_Groupmap}{$poplabel};
		push @{ $poplabels_bygroup{$grouplabel} }, $poplabel;
		push @{ $popindices_bygroup{$grouplabel} }, $self->{_Popidx}{$poplabel};
		}
	
	for my $grouplabel (@{ $self->{_Grouplabels} }){
	    if (exists $poplabels_bygroup{$grouplabel}){
			push @subgrouplabels, $grouplabel;
			push @popgroups, $popindices_bygroup{$grouplabel};
			}
		}

  	$self->changeGrouplabels(\@subgrouplabels);
	$self->{_Groups}=\@popgroups;
	$self->{_Groupsizes}=[ map { Misc::sum(@{ $self->{_Popsizes} }[@$_]) } @popgroups ];
	$self->{_Poplabels_bygroup}=\%poplabels_bygroup;
	return;
	}


sub changePoplabels{
	my ($self, $ref_poplabels, $ref_popsizes, $ref_sampleids_bypop)=@_;
	unless (@_==4){ die "ERROR: Wrong number of arguments!\n" }
	unless(@$ref_poplabels==@$ref_popsizes){ die "ERROR: Population label and population size array don't match!\n" }
	unless(@$ref_poplabels==@$ref_sampleids_bypop){ die "ERROR: Population label and individualIDs per pop array don't match!\n" }
  	$self->{_Poplabels}=$ref_poplabels;
	$self->{_Popidx}={ map { $ref_poplabels->[$_] => $_ } (0..$#{ $ref_poplabels }) };
	$self->{_Popsizes}=$ref_popsizes;
	$self->{_Sampleids_bypop}=$ref_sampleids_bypop;
	return;
	}

sub changeGrouplabels{
	my ($self, $ref_grouplabels)=@_;
	unless (@_==2){ die "ERROR: Wrong number of arguments!\n" }
  	$self->{_Grouplabels}=$ref_grouplabels;
	$self->{_Groupidx}={ map { $ref_grouplabels->[$_] => $_ } (0..$#{ $ref_grouplabels }) };
	return;
	}


1;

