#!/usr/bin/env perl

use strict;
use warnings;
use Classes::Populations;
use Classes::Fasta;
use Classes::Regions;
use Classes::Missing;
use Classes::SFS;
use Classes::Misc;

# written by Alexander Nater, March 2018
if (@ARGV < 13){
	die "Wrong number of arguments!\nUsage: ./generate_validbed_Pongo.pl [reference fasta file] [ancestral reference fasta file or 'none'] [bed file of regions to be excluded or 'none'] [output file] [allsiteslist file] [poplist file] [grouplist file or 'none'] [comma separated list of chromosome names or 'all' if working over whole genome] [comma separated list of population labels to include or 'all' if all populatons should be included] [comma separated list of individual ids or 'all'] [minimum depth per individual] [max missingness on population and group level] [mask regions within x bases of excluded regions, 0 for no masking] [verbose output 1=yes/0=no]\n";
	}
print STDERR $_+1, ". Argument: ", $ARGV[$_], "\n" foreach (0..$#ARGV);

# Set your options here:
my $maxsize=500000;	# combined maximum size of regions to process together.
my $minmappingquality=20;
my $minvariantquality=0;

my $fasta_file=shift @ARGV;
my $ancref_file=shift @ARGV;
my $excluded_file=shift @ARGV;
my $outfile=shift @ARGV;
my $allsiteslist=shift @ARGV;
my $poplist=shift @ARGV;
my $grouplist=shift @ARGV;
my @ids=Misc::parseArglist(shift @ARGV);
my @selpoplabels=split(',', shift @ARGV);
my @selindividuals=split(',', shift @ARGV);
my $mincoverage=shift @ARGV;
my $missingness=shift @ARGV;
my $excluded_dist=shift @ARGV;
my $verbose=shift @ARGV;

#-----------------------------------------------------------------------------------------------------------------------------

# initialize population object:
my $populations=Populations->new($poplist, $grouplist, \@selpoplabels, \@selindividuals);

# prepare reference fasta file:
my $fasta=Fasta->new();
$fasta->readIndex($fasta_file . ".fai", $fasta_file, 0);
my $ancref=undef;
if ($ancref_file ne "none"){
    $ancref=Fasta->new();
    $ancref->readIndex($ancref_file . ".fai", $ancref_file, 0);
    }

# prepare windows.
my $regions=$fasta->locifromIndex(undef, @ids);
my $windows=$regions->windowingLoci(0, $maxsize, $maxsize);

# exclude windows near genic regions:
if ($excluded_dist && $excluded_file ne "none"){
    my $excluded=Regions->new();
    $excluded->readBED($excluded_file, 'bed', 0);
    $windows=$windows->splitRegions($excluded, $excluded_dist, $verbose);
    $windows->printLoci("-");
    }

# prepare SFS object:
my $sfs=SFS->new();
my $validsites=0;

# initialize missing data object:
my $missing=Missing->new($allsiteslist, 'vcf');
$missing->setPopulations($populations, $verbose);

# prepare array of minimum number of individuals per population:
my $ref_subpoplabels=$populations->getPoplabels();
my $ref_nindividuals=$populations->getPopsizes();
my @minpopind=map { $_-int($missingness*$_) } @$ref_nindividuals;
print "Populations: ", join(',', @$ref_subpoplabels), "\n";
print STDERR "Number of individuals per population: ", join(',', @$ref_nindividuals), "\n";
print STDERR "Minimum number of individuals per population: ", join(',', @minpopind), "\n";

# prepare array of minimum number of individuals per group:
my $ref_popgroups=$populations->getGroups();
my $ref_groupsizes=$populations->getGroupsizes();
my @mingroupind=map { $_-int($missingness*$_) } @$ref_groupsizes;
print STDERR $populations->getGrouplabels($_), ": ", join(',', @{ $ref_popgroups->[$_] }), "\n" foreach (0..$#{ $ref_popgroups });
print STDERR "Number of individuals per group: ", join(',', @$ref_groupsizes), "\n";
print STDERR "Minimum number of individuals per group: ", join(',', @mingroupind), "\n";


SUBSET: while (1){
	# get a subset of loci with a combined size smaller than maxsize:
	my $subset=$windows->subsetLoci($maxsize);
	unless ( $subset->getSize() ){ last SUBSET }
	$subset->printLoci("-");

	# acquire missing data:
	$missing->readAllsites($subset, $minmappingquality, $minvariantquality, undef, 0, 0, $verbose);

	# generate bed file of valid regions:
	$validsites+=$sfs->validRegions($outfile, $subset, $fasta, $ancref, $missing, $ref_nindividuals, \@minpopind, $ref_popgroups, \@mingroupind, $mincoverage, $verbose);
	}

print STDERR "Found a total of $validsites valid sites.\n";

