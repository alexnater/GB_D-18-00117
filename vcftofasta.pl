#!/usr/bin/env perl

use strict;
use warnings;
use Classes::Populations;
use Classes::Fasta;
use Classes::Regions;
use Classes::Missing;
use Classes::Genotypes;
use Classes::Misc;

# written by Alexander Nater, July 2016
if (@ARGV < 15){
    die "Wrong number of arguments!\nUsage: ./vcftofasta_Pongo.pl [reference fasta file] [output folder] [outfile prefix] [allsiteslist file] [vcflist file] [poplist file] [grouplist file or 'none'] [comma separated list of chromosome names or 'all' if working over whole genome] [comma separated list of population labels to include or 'all' if all populatons should be included] [comma separated list of individuals or 'all' to be printed] [bed file with predefined regions or 'none' for whole genome] [minimum depth per individual] [set reference Ns to missing?] [haplotize individuals?] [max missingness on population and group level] [verbose output 1=yes/0=no]\n";
    }
print STDERR $_+1, ". Argument: ", $ARGV[$_], "\n" foreach (0..$#ARGV);


# Set your options here:
my $maxsize=5000000;    # combined maximum size of regions to process together.
my $minmappingquality=20;
my $minvariantquality=0;
my $mingenotypequality=0;    # this setting biases towards lower diversity if > 0!
my $rowlength=80;

my $fasta_file=shift @ARGV;
my $outfolder=shift @ARGV;
my $outprefix=shift @ARGV;
my $allsiteslist=shift @ARGV;
my $vcflist=shift @ARGV;
my $poplist=shift @ARGV;
my $grouplist=shift @ARGV;
my @ids=Misc::parseArglist(shift @ARGV);
my @selpoplabels=split(',', shift @ARGV);
my @selindividuals=split(',', shift @ARGV);
my $locifile=shift @ARGV;
my $mincoverage=shift @ARGV;
my $excluderefN=shift @ARGV;    # FALSE=0/TRUE=1
my $haplotize=shift @ARGV;
my $missingness=shift @ARGV;
my $verbose=shift @ARGV;

#-----------------------------------------------------------------------------------------------------------------------------

$outfolder=~s{/\z}{};    # remove trailing slash from folder path.
mkdir $outfolder unless (-d "$outfolder");

# initialize population object:
my $populations=Populations->new($poplist, $grouplist, \@selpoplabels, ['all']);

# prepare reference fasta file:
my $fasta=Fasta->new();
$fasta->readIndex($fasta_file . ".fai", $fasta_file);

# get regions from bed file:
my $regions;
if (defined $locifile && $locifile ne 'none'){
    print "Reading loci from $locifile ...\n";
    $regions=Regions->new();
    $regions->readBED($locifile, 'bed', 0);
    }

# or get regions from fasta index:
else { $regions=$fasta->locifromIndex() }
my $windows=$regions->windowingLoci(0, $maxsize, $maxsize, 0, @ids);

# initialize missing data object:
my $missing=Missing->new($allsiteslist, 'vcf');
$missing->setPopulations($populations, $verbose);

# initialize genotype data object:
my $genotypes=Genotypes->new($vcflist);
$genotypes->setPopulations($populations, $verbose);

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

# obtain array of individual indices to be printed:
my %subindlabels;
%subindlabels=map { $_=>1 } @selindividuals if ($selindividuals[0] ne 'all');
my (@subindividuals, @indnames);
my $totsize=0;
my $ref_indnames_bypop=$missing->getIndnames();
for my $ref_inds (@$ref_indnames_bypop){
    my @tempind;
    if ($selindividuals[0] eq 'all'){ @tempind=( 0..$#{ $ref_inds } ) }
    else {
        for my $index (0..$#{ $ref_inds }){
            push @tempind, $index if ($subindlabels{ $ref_inds->[$index] });    # test if list of subpoplabels contain individual ids.
            }
        }
    push @subindividuals, \@tempind;
    push @indnames, map { $ref_inds->[$_] } @tempind;
    $totsize+=@tempind;
    }
print STDERR join(',', @indnames), "\n";
die "Wrong number of individual names!\n" unless (@indnames==$totsize);

# prepare fasta object for each individual to be printed:
my @fasta_objects;
push @fasta_objects, Fasta->new() foreach (0..$totsize-1);


# process windows:
SUBSET: while (1){
    # get a subset of loci with a combined size smaller than maxsize:
    my $subset=$windows->subsetLoci($maxsize);
    unless ( $subset->getSize() ){ last SUBSET }
    $subset->printLoci("-");

    # acquire missing data:
    $missing->readAllsites($subset, $minmappingquality, $minvariantquality, undef, 0, 0, $verbose);

    # acquire genotype data:
    $genotypes->readvcfList($subset, $mingenotypequality, 0, 1, 0, 0, 0, $verbose);

    # obtain individual fasta sequences:
    $subset->generateFasta($fasta, \@fasta_objects, $genotypes, $missing, \@minpopind, $ref_popgroups, \@mingroupind, $mincoverage, $excluderefN, $haplotize, \@subindividuals);
    }


# print fasta sequence for each selected individual:
for my $ind (0..$totsize-1){
    my $indfolder=$indnames[$ind];
    mkdir "/$outfolder/$indfolder" unless (-d "/$outfolder/$indfolder");
    my $outname="$outfolder" . "/" . "$indfolder" . "/" . "$outprefix" . "_" . "$indnames[$ind]" . "_" . $ids[0] . "-" . $ids[-1] . ".fasta";
    $fasta_objects[$ind]->printSeqs($outname, $rowlength);
    }


