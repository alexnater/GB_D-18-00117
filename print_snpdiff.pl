#!/usr/bin/env perl

use strict;
use warnings;
use Classes::Populations;
use Classes::Fasta;
use Classes::Regions;
use Classes::Genotypes;
use Classes::Coding;
use Classes::SFS;

# written by Alexander Nater, July 2016
if (@ARGV < 15){
    die "Wrong number of arguments!\nUsage: ./pring_snpdiff_Pongo.pl [reference fasta file] [ancestral fasta file] [annotation gtf file] [output folder] [vcflist file] [poplist file] [grouplist file or 'none'] [comma separated list of chromosome names or 'all' if working over whole genome] [comma separated list of population labels to include or 'all' if all populatons should be included] [comma separated list of individuals or 'all'] [bed file with predefined regions or 'none' for whole genome] [minimum depth per individual] [print only annotated SNPs?] [haplotize individuals?] [max missingness on population and group level] [verbose output 1=yes/0=no]\n";
    }
print STDERR $_+1, ". Argument: ", $ARGV[$_], "\n" foreach (0..$#ARGV);


# Set your options here:
my $mingenotypequality=0;    # this setting biases towards lower diversity if > 0!
my $maxsize=10000000;

my $fasta_file=shift @ARGV;
my $ancref_file=shift @ARGV;
my $genes_file=shift @ARGV;
my $outfolder=shift @ARGV;
my $vcflist=shift @ARGV;
my $poplist=shift @ARGV;
my $grouplist=shift @ARGV;
my @ids=Misc::parseArglist(shift @ARGV);
my @selpoplabels=split(',', shift @ARGV);
my @selindividuals=split(',', shift @ARGV);
my $locifile=shift @ARGV;
my $mincoverage=shift @ARGV;
my $onlyannotated=shift @ARGV;    #FALSE=0/TRUE=1
my $haplotize=shift @ARGV;    #FALSE=0/TRUE=1
my $missingness=shift @ARGV;
my $verbose=shift @ARGV;

#-----------------------------------------------------------------------------------------------------------------------------

$outfolder=~s{/\z}{};    # remove trailing slash from folder path.
mkdir $outfolder unless (-d "$outfolder");
my $idstring=@ids>1 ? $ids[0] . "_" . $ids[-1] : $ids[0];
my $outfile="diffreport_" . $idstring . ".txt";
my $annotation_report="annotation_report_" . $idstring . ".txt";

# initialize population object:
my $populations=Populations->new($poplist, $grouplist, \@selpoplabels, \@selindividuals);

# prepare reference fasta file:
my $fasta=Fasta->new();
$fasta->readIndex($fasta_file . ".fai", $fasta_file);
my $ancref=Fasta->new();
$ancref->readIndex($ancref_file . ".fai", $ancref_file, 0);

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

# prepare locus definitions of coding regions:
my $cds=Coding->new();
$cds->readCDS($genes_file, 0, undef, \@ids);
my $genes=$cds->concatenateTranscripts($fasta, $ancref, 0, 2);
if (1){     # restricts live of $exon_genotypes object to block.
    my $concat_genes=$genes->concatRegions();
    $genotypes->readvcfList($concat_genes, $mingenotypequality, 0, 1, 0, 0, 0, $verbose);
    $cds->annotateSNPs($genes, $genotypes, "$outfolder/$annotation_report");
    }
$genes->printLoci("-");

# create SFS object:
my $sfs=SFS->new();


SUBSET: while (1){
    # get a subset of loci with a combined size smaller than maxsize:
    my $subset=$windows->subsetLoci($maxsize);
    $subset->printLoci("-");
    unless ( $subset->getSize() ){ last SUBSET }

    # acquire genotype data:
    $genotypes->readvcfList($subset, $mingenotypequality, 0, 1, 1, 0, 0, $verbose);

    # print report of frequency differences:
    $sfs->printSNPDifferentials("$outfolder/$outfile", $subset, $cds, $fasta, $ancref, undef, $genotypes, $ref_nindividuals, $ref_popgroups, \@mingroupind, $mincoverage, $haplotize, undef, $onlyannotated);
    }


