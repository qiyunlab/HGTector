#!/usr/bin/perl

use warnings;
use strict;
$| = 1;

print "
-> Analyzer: Identify putative HGT-derived genes based on search results. <-
";

print "

This program analyzes search results of proteins from each genome
and identify ones that may have undergone horizontal gene transfer
events and their donor and/or recipients, as well as any potential
gene origination and loss events.

Every subfolder in search/ is considerd as protein set. The result
is saved as a tab-delimited file, which can be visualized by programs
such as Excel.

Usage:
  perl analyzer.pl <working directory>

Output:
  result/<name>.txt

" and exit unless @ARGV;



## all-purpose variables ##

my $i; my $j; my $n; my $s; my $t; my @a; my @b; my @c; my %h;


## subroutines ##

sub median(@);
sub mad(@);
sub quantiles(@);
sub z_scores(@);
sub z_test(@);
sub modified_z(@);
sub boxplot(@);
sub adjusted_boxplot(@);
sub recurse_deOutlier(@);
sub recurse_Z(@);


## global variables ##

my @sets;                                        # proteins sets to analyze
my %lv = ();                                            # taxonomic grouping scenario
my %selves = ();                                # taxonomic information of self group
my %lvList = ();                                # list of taxids of each level

my %taxadb = ();                                # taxa.db
my %ranksdb = ();                                # ranks.db
my %selfinfo = ();                                # self.info

my %proteins = ();                                # accn -> name (identified by COG)

my @files;                                        # search report files for each protein set

## the master table storing everything. It's an array of hashes. Each row is a record.
my %results = ();

## the phyletic pattern of the whole genome, aka "fingerprint"
my %fpN = ();                                    # number of hits per protein

## program parameters ##

my $wkDir = $ARGV[0];                            # working directory
my $interactive = 1;                            # interactive or automatic mode

my $minHits = 0;                                # minimal number of hits a valid search result should contain
my $maxHits = 0;                                # maximal number of hits to retain for one protein, 0 means infinite
my $minSize = 0;                                # minimal size (aa) of a valid protein (0 means infinite)
my $evalue = 1e-5;                                # E-value cutoff
my $identity = 0;                                # percent identity cutoff
my $coverage = 0;                                # query coverage cutoff

# algorithm
my $selfRank = 0;                                # taxonomic rank(s) on which the program analyzes
my $normalize = 1;                                # use relative bit score (bit score of subject / bit score of query)
my $unite = 1;                                    # hit pattern (0: each genome has own pattern, 1: one pattern for all genomes)

my $useDistance = 0;                            # use phylogenetic distance instead of bit scores
my $useWeight = 1;                                # use weight (sum of scores) instead of number of hits

# fingerprints
my $outRaw = 1;                                    # output raw number/weight data
my $outFp = 1;                                    # output fingerprint
my $graphFp = 0;                                # graph fingerprint (requires R)

my $plotRef = "";                                # file name of reference set of genes

my $boxPlot = 1;                                # box plot
my $histogram = 1;                                # histogram
my $densityPlot = 1;                            # density plot
my $scatterPlot = 1;                            # scatter plot
my $plot3D = 0;                                    # 3-way scatter plot

# cutoffs
my $howCO = 4;                                    # how to determine cutoffs (0: user-defined global cutoff (%), 1: user-defined individual
                                                    # cutoffs, 2: wait for user input, 3: histogram, 4: kernel density estimation,
                                                    # 5: hierarchical clustering)
my $globalCO = 0.25;                            # arbitrary global cutoff (%)
my ($selfCO, $closeCO, $distalCO) = (0, 0, 0);    # user-defined cutoffs for individual groups

my $exOutlier = 0;                                # exclude outliers, hits distant from other hits of the same group

my $nBin = 20;                                    # number of bins in histogram
my $plotHist = 0;                                # plot histogram on screen

my $toolKDE = 0;                                # computational tool for kernel density estimation
my $bwF = 1;                                    # bandwidth selection factor
my $plotKDE = 0;                                # plot density function on screen
my $toolExtrema = 0;                            # computational tool for identifying local extrema of density function (0: Perl code, 1: R package "pastecs")
my $whichPeak = 0;                                # definition of "typical" and "atypical" regions
my $modKCO = 1;                                    # location of cutoff (0: 1st pit, 1: midpoint of x-coordinates between 1st peak and 1st pit, 2/3: quantile
my $qKCO = 0.5;                                    # horizontal/vertical quantile from 1st pit toward 1st peak

my $dipTest = 0;                                # perform non-unimodality test (Hartigan's dip test) and report p-value
my $dipSig = 0;                                    # use global cutoff if dip test's result is not significant

my $selfLow = 0;                                # HGT-derived genes must have low self weight (an optional criterion)

my $BBH = 0;                                    # use conventional best match method instead
my $loss = 0;                                    # also report gene loss events
my $POE = 0;                                    # also report POE

my @ranks = ('species', 'genus', 'family', 'order', 'class', 'phylum');

my @selfGroup = ();
my @closeGroup = ();
my @excludeGroup = ();
my @inSets = ();
my @exSets = ();

my $R;                                            # Statistics::R instance

## read configurations ##

if (-e "$wkDir/config.txt") {
    open IN, "<$wkDir/config.txt";
    while (<IN>) {
        s/#.*$//; s/\s+$//g; s/^\s+//g; next unless $_;
        $interactive = $1 if /^interactive=([01])$/;
        $minHits = $1 if /^minHits=(\d+)$/;
        $maxHits = $1 if /^maxHits=(\d+)$/;
        $minSize = $1 if /^minSize=(\d+)$/;
        $evalue = $1 if /^evalue=(.+)$/;
        $identity = $1 if /^identity=(.+)$/;
        $identity = $1 if /^percIdent=(.+)$/; # backward compatibility
        $coverage = $1 if /^coverage=(.+)$/;

        $selfRank = $1 if /^selfRank=(.+)$/;
        $normalize = $1 if /^normalize=(.+)$/;
        $unite = $1 if /^unite=([01])$/;
        $useDistance = $1 if /^useDistance=([01])$/;
        $useWeight = $1 if /^useWeight=([01])$/;

        $howCO = $1 if /^howCO=(\d)$/;
        $globalCO = $1 if /^globalCO=(.+)$/;
        $selfCO = $1 if /^selfCO=(.+)$/;
        $closeCO = $1 if /^closeCO=(.+)$/;
        $distalCO = $1 if /^distalCO=(.+)$/;

        $exOutlier = $1 if /^exOutlier=([0123])$/;

        $nBin = $1 if /^nBin=(\d+)$/;
        $plotHist = $1 if /^plotHist=([01])$/;

        $toolKDE = $1 if /^toolKDE=([01])$/;
        $bwF = $1 if /^bwF=(.+)$/;
        $plotKDE = $1 if /^plotKDE=([01])$/;
        $toolExtrema = $1 if /^toolExtrema=([01])$/;
        $whichPeak = $1 if /^whichPeak=([0123])$/;
        $modKCO = $1 if /^modKCO=([0123])$/;
        $qKCO = $1 if /^qKCO=(.+)$/;

        $dipTest = $1 if /^dipTest=([01])$/;
        $dipSig = $1 if /^dipSig=(.+)$/;

        $selfLow = $1 if /^selfLow=([01])$/;

        $outRaw = $1 if /^outRaw=([01])$/;
        $outFp = $1 if /^outFp=([01])$/;
        $graphFp = $1 if /^graphFp=([01])$/;
        $plotRef = $1 if /^plotRef=(.+)$/;

        $boxPlot = $1 if /^boxPlot=([01])$/;
        $histogram = $1 if /^histogram=([01])$/;
        $densityPlot = $1 if /^densityPlot=([01])$/;
        $scatterPlot = $1 if /^scatterPlot=([01])$/;
        $plot3D = $1 if /^plot3D=([01])$/;

        $BBH = $1 if /^BBH=([012])$/;
        $loss = $1 if /^loss=([01])$/;
        $POE = $1 if /^POE=([01])$/;

        @ranks = split (/\s*,\s*/, $1) if /^ranks=(.+)$/;

        @selfGroup = split (/\s*,\s*/, $1) if /^selfGroup=(.+)$/;
        @closeGroup = split (/\s*,\s*/, $1) if /^closeGroup=(.+)$/;
        @excludeGroup = split (/\s*,\s*/, $1) if /^excludeGroup=(.+)$/;
        @inSets = split (/\s*,\s*/, $1) if /^inSets=(.+)$/;
        @exSets = split (/\s*,\s*/, $1) if /^exSets=(.+)$/;
    }
    close IN;
}
if ($identity and $identity < 1) { $identity *= 100; }
if ($coverage and $coverage < 1) { $coverage *= 100; }

## check previous result ##

if (-d "$wkDir/result/detail" and -d "$wkDir/result/statistics") {
    print "Warning: Prediction result from a previous analysis is detected.\n";
    if ($interactive) { print "Press Enter to overwrite, or Ctrl+C to exit:\n"; $s = <STDIN>; }
    else { print "To be overwritten.\n"; }
}

## verify configurations ##

if ($globalCO) {
    $globalCO = $globalCO / 100 if ($globalCO =~ s/%$//);
    die "Error: Global cutoff must be between 0 and 1.\n" if ($globalCO <= 0 or $globalCO >=1);
}

# conditionally use Statistics::R for Perl-R communication;
if ($graphFp or ($howCO == 4 and ($toolKDE or $toolExtrema)) or ($howCO == 5) or $dipTest) {
    eval { require Statistics::R; Statistics::R->import() };
    die "Error: Perl module Statistics::R is not available.\n" if ($@);
    $R = Statistics::R->new();
}

## initiate global variables ##

if ($unite) {
    @b = (); @a = ('0','1','2');
    $fpN{'0'}{$_}{'data'} = [@b] for (@a);
    # $fpS{'0'}{$_}{'data'} = [@b] for (@a);
}


## Identify taxonomic levels ##

print "Reading taxonomic information...";
open IN, "<$wkDir/taxonomy/taxa.db";
while (<IN>) {
    next if /^#/;
    @a = split /\t/;
    next unless $#a;
    $a[$#a] =~ s/\s+$//;
    %h = ('name',$a[1],'rank',$a[2]);
    $i = 3; $h{$_} = $a[$i++] for (@ranks);
    $taxadb{$a[0]} = {%h};
}
close IN;
open IN, "<$wkDir/taxonomy/ranks.db";
while (<IN>) {
    s/\s+$//; next if /^#/; next unless $_;
    @a = split /\t/;
    $ranksdb{$a[0]} = $a[1];
}
close IN;
open IN, "<$wkDir/taxonomy/self.info";
while (<IN>) {
    s/\s+$//; next if /^#/; next unless $_;
    @a = split /\t/;
    next if (scalar(@a) < 3);
    %h = ('taxid',$a[1],'name',$a[2]);
    $selfinfo{$a[0]} = {%h};
}
close IN;
print " done.\n";

print "Analyzing taxonomic information...";
for ($i=0; $i<=$#ranks; $i++) {
    $j = 1; # if all names are the same
    $s = 0; # last name of this rank
    foreach my $key (keys %selfinfo) {
        $t = $taxadb{$selfinfo{$key}{'taxid'}}{$ranks[$i]};
        $s = $t unless $s;
        if ($t != $s) { $j = 0; last; }
    }
    last if $j;
}
print " done.\n";
$selves{'rank'} = $ranks[$i];
$selves{'rank_id'} = $i;
$selves{'rank_taxid'} = $s;
if (exists $ranksdb{$s}) { $selves{'rank_name'} = $ranksdb{$s}; }
elsif (exists $taxadb{$s}) { $selves{'rank_name'} = $taxadb{$s}{'name'}; }
else { die "Unknown TaxID $s.\n"; }
print "  All input genomes belong to $selves{'rank'} $selves{'rank_name'} (TaxID: $selves{'rank_taxid'}).\n";

if (@selfGroup) { # user-defined self group
    $lv{'rank'} = "(user-defined self)"; $lv{'id'} = -1; $lv{'taxid'} = join (",", @selfGroup);
    @a = ();
    foreach (@selfGroup) {
        if (exists $ranksdb{$_}) { push @a, $ranksdb{$_}; }
        elsif (exists $taxadb{$_}) { push @a, $taxadb{$_}{'name'}; }
        else { push @a, "unknown"; }
    }
    $lv{'name'} = join (",", @a);
} else { # lowest rank that contains all input organisms
    if ($selfRank =~ /\d/) { # relative level number
        $lv{'rank'} = $ranks[$selves{'rank_id'} + $selfRank];
        $lv{'id'} = $selves{'rank_id'} + $selfRank;
    } else { # rank name
        $lv{'rank'} = $selfRank;
        for ($j=0; $j<=$#ranks; $j++) {
            if ($ranks[$j] eq $selfRank) { $lv{'id'} = $j; last; }
        }
    }
    $s = $selfinfo{(keys %selfinfo)[0]}{'taxid'};
    $lv{'taxid'} = $taxadb{$s}{$lv{'rank'}};
    if (exists $ranksdb{$lv{'taxid'}}) { $lv{'name'} = $ranksdb{$lv{'taxid'}}; }
    elsif (exists $taxadb{$lv{'taxid'}}) { $lv{'name'} = $taxadb{$lv{'taxid'}}{'name'}; }
    else { die "Unknown TaxID $lv{'taxid'}.\n"; }
}
%h = ();
foreach my $key (keys %taxadb) {
    next unless $taxadb{$key}{'rank'};
    foreach (split (/,/, $lv{'taxid'})) {
        $h{$key} = 1 if ($key eq $_ or $taxadb{$key}{'rank'} =~ /\/$_$/ or $taxadb{$key}{'rank'} =~ /\/$_\//);
    }
}
$lvList{$lv{'rank'}} = {%h};

if (@closeGroup) { # user-defined close group
    $lv{'prank'} = "(user-defined close)"; $lv{'ptaxid'} = join (",", @closeGroup);
    @a = ();
    foreach (@closeGroup) {
        if (exists $ranksdb{$_}) { push @a, $ranksdb{$_}; }
        elsif (exists $taxadb{$_}) { push @a, $taxadb{$_}{'name'}; }
        else { push @a, "unknown"; }
    }
    $lv{'pname'} = join (",", @a);
} else { # lowest parent rank that has adequate members for statistical analysis
    if ($lv{'id'} >= $#ranks) { die "Cannot find a parent taxonomic rank of $lv{'taxid'}.\n"; }
    print "  Choose one of the following parental taxonomic ranks as the close group:\n";
    my @pranks = ();
    $s = $selfinfo{(keys %selfinfo)[0]}{'taxid'};
    for ($i=$lv{'id'}+1; $i<=$#ranks; $i++) {
        next unless exists $taxadb{$s}{$ranks[$i]};
        $t = $taxadb{$s}{$ranks[$i]};
        %h = ('rank' => $ranks[$i], 'taxid' => $t, 'name' => '', 'n' => 0, 'm' => 0); # m is the number of members - number of self members
        if (exists $ranksdb{$t}) { $h{'name'} = $ranksdb{$t}; }
        elsif (exists $taxadb{$t}) { $h{'name'} = $taxadb{$t}{'name'}; }
        else { next; }
        foreach my $key (keys %taxadb) {
            if ($key eq $t or $taxadb{$key}{'rank'} =~ /\/$t$/ or $taxadb{$key}{'rank'} =~ /\/$t\//) {
                $h{'n'} ++;
                $h{'m'} ++ unless exists $lvList{$lv{'rank'}}{$key};
            }
        }
        push @pranks, {%h};
        print "    ".$h{'rank'}." ".$h{'name'}." (TaxID: ".$h{'taxid'}.") (".$h{'n'}." members).\n";
    }
    unless (scalar @pranks) { die "Cannot find a parent taxonomic rank of $lv{'taxid'}.\n"; }
    for ($i=0; $i<=$#pranks; $i++) {
        if ($pranks[$i]{'m'} >= 10 or $i == $#pranks) { ######## new rule: close group should have >=10 members
            print "  The program intelligently chose $pranks[$i]{'rank'} $pranks[$i]{'name'}.\n";
            if ($interactive) { print "Press Enter to accept, or Ctrl+C to exit:\n"; $s = <STDIN>; }
            $lv{'prank'} = $pranks[$i]{'rank'};
            $lv{'ptaxid'} = $pranks[$i]{'taxid'};
            $lv{'pname'} = $pranks[$i]{'name'};
            last;
        }
    }
}
%h = ();
foreach my $key (keys %taxadb) {
    next unless $taxadb{$key}{'rank'};
    foreach (split (/,/, $lv{'ptaxid'})) {
        $h{$key} = 1 if ($key eq $_ or $taxadb{$key}{'rank'} =~ /\/$_$/ or $taxadb{$key}{'rank'} =~ /\/$_\//);
    }
}
$lvList{$lv{'prank'}} = {%h};

print "Analysis will work on the following taxonomic ranks:\n";
print "  Self: $lv{'rank'} $lv{'name'} (TaxID: $lv{'taxid'}) (".(keys %{$lvList{$lv{'rank'}}})." members),\n";
print "  Close: $lv{'prank'} $lv{'pname'} (TaxID: $lv{'ptaxid'}) (".(keys %{$lvList{$lv{'prank'}}})." members),\n";
print "  Distal: all other organisms.\n";
if ($interactive) {
    print "Press Enter to continue, or Ctrl+C to exit:";
    $s = <STDIN>;
}

## Read reference proteins ##

my %refs = ();
if ($plotRef) {
    if (-e $plotRef) { open IN, "<$plotRef"; }
    elsif (-e "$wkDir/$plotRef") { open IN, "<$wkDir/$plotRef"; }
    else { die "Reference file $plotRef does not exist.\n"; }
    while (<IN>) {
        s/\s+$//; next unless $_;
        @a = split (/\t/);
        $a[0] =~ s/\.\d+$//;
        $refs{$a[0]} = [0,0,0,0]; # exist, self, close, distal
    }
    close IN;
}


## Read protein sets ##

print "Reading protein sets...";
opendir (DIR, "$wkDir/search");
@a = readdir(DIR);
close DIR;
foreach (@a) {
    next if (/^\./);
    push @sets, $_ if -d "$wkDir/search/$_";
}
print "done. ";
die "No genome detected.\n" unless @sets;
print @sets." sets detected.\n";


# Summarize search reports ##

print "Analyzing search results...\n";
print "0-------------25-------------50------------75------------100%\n";

foreach my $set (@sets) {
    if (@inSets) { $i = 0; foreach (@inSets) { if ($set eq $_) { $i = 1; last; } } next unless $i; }
    if (@exSets) { $i = 0; foreach (@exSets) { if ($set eq $_) { $i = 1; last; } } next if $i; }
    opendir (DIR, "$wkDir/search/$set");
    @files = grep(/\.txt$/,readdir(DIR));
    close DIR;
    print "No protein found in $set\n" and next unless @files;

    ## varibles to show a progress bar
    my $iProtein = 0;
    my $iProgress = 0;
    my $nProtein = $#files+1;
    print "$set has $nProtein proteins. Analyzing...\n";

    unless ($unite) {
        @a = ('0','1','2'); @b = ();
        $fpN{$set}{$_}{'data'} = [@b] for (@a);
    }

    ## information of self
    my %self = ();

    foreach my $file (@files) {
        $iProtein ++;

        my %result = ();                # a record to store everything about this search
        my %scores = ();                # scores of each hit by category, as a buffer for computing the statistics above
        my @hits;                        # parameters of the hits. one hit contains:
                                        # accn, organism, group, taxid, genus, score
        $file =~ /(.+)\.[^.]+$/;
        $result{'query'} = $1;

        my $nHits = 0;                    # total number of hits. just for convenience
        my $nScore = 0;                    # total score
        my $lastHit = "";                # store the last hit in the organism table, for the identification of duplicated taxa
        my $selfScore = 0;

        # read hit table #

        open IN, "<$wkDir/search/$set/$file" or next;
        my $reading = 0;
        my ($hasCoverage, $hasDistance) = (0, 0);
        while (<IN>) {
            s/\s+$//;
            if (/^BEGIN QUERY/) { $reading = "query"; next; }
            if (/^BEGIN ORGANISM/) { $reading = "organism"; next; }
            if (/^BEGIN DATA/) { $reading = "data"; next; }
            if (/^END;/) { $reading = 0; next; }
            if ($reading eq "query") { # read query (self)
                $result{'accn'} = $1 if /^\tName=(.+);$/;
                $result{'gi'} = $1 if /^\tGI=(\d+);$/;
                $result{'length'} = $1 if /^\tLength=(.+);$/;
                $result{'product'} = $1 if /^\tProduct=(.+)\s*;$/;
                $result{'organism'} = $1 if /^\tOrganism=(.+)\s*;$/;
                if (/^\tAccession=(.+);$/) { $result{'accn'} = $1; $result{'accn'} =~ s/\.[\d]+$//; }
            }
            if ($reading eq "organism") { # read organisms
                next if /^;/;
                if (/^\[/) {
                    $hasCoverage = 1 if /Coverage/;
                    $hasDistance = 1 if /Distance/;
                    next;
                }
                @a = split (/\t/);
                if ($#a < 5) {
                    print "\nIncomplete hit record $a[0] in $file of $set.\n";
                    if ($interactive) {
                        print "Press Enter to continue, or Ctrl+C to exit:";
                        $s = <STDIN>;
                    }
                    next;
                }

                # filter out low-quality hits
                next if ($a[$#a] eq "x");
                next if ($evalue and $a[4] ne "*" and $a[4] > $evalue);
                next if ($identity and $a[5] ne "*" and $a[5] < $identity);
                next if ($coverage and $hasCoverage and $a[6] and $a[6] ne "*" and $a[6] < $coverage);

                # filter out unwanted taxonomy groups
                if (@excludeGroup) {
                    my $isExclude = 0;
                    foreach (@excludeGroup) {
                        if ($a[2] == $_ or $taxadb{$a[2]}{'rank'} =~ /\/$_$/ or $taxadb{$a[2]}{'rank'} =~ /\/$_\//) {
                            $isExclude = 1;
                            last;
                        }
                    }
                    next if $isExclude;
                }

                my %hit = ();
                $hit{'accns'} = $a[0];
                $hit{'organism'} = $a[1];
                $hit{'taxid'} = $a[2];
                $hit{'score'} = $a[3];
                $hit{'evalue'} = $a[4];
                $hit{'identity'} = $a[5];
                $hit{'coverage'} = $a[6] if $hasCoverage;

                # use phylogenetic distance instead of bit score
                if ($useDistance) {
                    if ($#a >= 6) { $hit{'score'} = 1 - $a[6+$hasCoverage]; }
                    elsif (!@hits) { $hit{'score'} = 0; }
                    else { $hit{'score'} = $hits[$#hits]{'score'}; }
                }

                @a = split(/\//, $a[0]);
                $hit{'accn'} = $a[0];
                push @hits, {%hit};
            }
            last if ($maxHits and $#hits >= $maxHits-1);
        }
        close IN;

        # skip if there is no hit
        next unless @hits;

        # next if ($minSize and ($result{'length'} < $minSize));
        unless (exists $result{'query'} and exists $result{'length'}) {
            print "\nIncomplete search result: $set/$file.\n" ;
            if ($interactive) {
                print "Press Enter to continue, or Ctrl+C to exit:";
                $s = <STDIN>;
            }
        }
        $result{'product'} = '' unless exists $result{'product'};

        ## Intepret hit table ##

        # total number of hits
        $result{'n'} = $#hits+1;

        # sort hits by bit score or phylogenetic distance
        @hits = sort {$b->{'score'} <=> $a->{'score'}} @hits;

        # identify self (query) information
        my $isQuery = 0;
        for ($i=0; $i<=$#hits; $i++) {
            @a = split(/\//, $hits[$i]{'accns'});
            foreach (@a) {
                if ($result{'accn'} eq $_) {
                    $result{'id'} = $i;
                    $result{'taxid'} = $hits[$i]{'taxid'};
                    $result{'score'} = $hits[$i]{'score'};
                    $result{'organism'} = $hits[$i]{'organism'};
                    $isQuery = 1;
                    last;
                }
            }
            last if exists $result{'id'};
        }
        unless (exists $result{'id'}) {
            $result{'id'} = 0;
            $result{'taxid'} = $hits[0]{'taxid'};
            $result{'score'} = $hits[0]{'score'};
            $result{'organism'} = $hits[0]{'organism'};
        }
        next unless $result{'score'};

        # Use absolute or relative bit scores

        if ($normalize and not $useDistance) {
            for ($i=0; $i<=$#hits; $i++) {
                $hits[$i]{'score'} = sprintf("%.3f", $hits[$i]{'score'}/$result{'score'});
            }
        }

        # initialize values of prediction results
        $result{'in'} = "";                    # whether incoming HGT or origination took place within the group
        $result{'loss'} = "";                # gene loss event
        $result{'origin'} = "";                # gene origination event
        $result{'income'} = "";                # incoming HGT event
        $result{'outcome'} = "";            # outcoming HGT event

        # Summarize numbers and scores #
        ## 0 - self group, 1 - close groups, 2 - distal gs
        ## N - number, S - scores
        ## hit1 - first close hit, hit 2 - first distal hit

        my ($topCloseScore, $topDistalScore) = (0, 0);
        $result{'N0'} = 0; $result{'N1'} = 0; $result{'N2'} = 0;
        @a = (); $result{'S0'} = [@a]; $result{'S1'} = [@a]; $result{'S2'} = [@a];
        for ($i=0; $i<=$#hits; $i++) {
            if (exists $lvList{$lv{'rank'}}{$hits[$i]{'taxid'}}) {
                next unless $isQuery;
                if ($useWeight) { $result{'N0'} += $hits[$i]{'score'}; }
                else { $result{'N0'} ++; }
                push @{$result{'S0'}}, $hits[$i]{'score'};
            } elsif (exists $lvList{$lv{'prank'}}{$hits[$i]{'taxid'}}) {
                if ($useWeight) { $result{'N1'} += $hits[$i]{'score'}; }
                else { $result{'N1'} ++; }
                push @{$result{'S1'}}, $hits[$i]{'score'};
                $result{'hit1'} = $hits[$i]{'taxid'} unless exists $result{'hit1'};
                $topCloseScore = $hits[$i]{'score'} unless $topCloseScore;
                $result{'BBH'} = "" unless exists $result{'BBH'};
            } else {
                if ($useWeight) { $result{'N2'} += $hits[$i]{'score'}; }
                else { $result{'N2'} ++; }
                push @{$result{'S2'}}, $hits[$i]{'score'};
                unless (exists $result{'hit2'}) {
                    $result{'hit2'} = $hits[$i]{'taxid'};
                    $result{'hit2accn'} = $hits[$i]{'accn'};
                }
                $topDistalScore = $hits[$i]{'score'} unless $topDistalScore;
                $result{'BBH'} = 1 unless exists $result{'BBH'};
            }
        }
        # more stringent BBH method:
        # $result{'BBH'} = 1 if ($topDistalScore > $topCloseScore);

        $result{'N0'} = sprintf("%.3f", $result{'N0'});
        $result{'N1'} = sprintf("%.3f", $result{'N1'});
        $result{'N2'} = sprintf("%.3f", $result{'N2'});

        if (exists $refs{$result{'accn'}}) {
            $refs{$result{'accn'}} = [1,$result{'N0'},$result{'N1'},$result{'N2'}];
        }

        # proteins without non-self hits are considered as de novo originated and not considered for statistics
        $result{'origin'} = 1 if $result{'N1'}+$result{'N2'} < 0.000001;

        # proteins with hits lower than minimum cutoff are considered as de novo originated and not considered for statistics
        if ($minHits and ($result{'n'} < $minHits)) { $result{'origin'} = 1; $result{'BBH'} = ""; }

        # Put this record into the master record
        push @{$results{$set}}, {%result};

        # Record fingerprint #
        if ($unite) { $s = '0'; }
        else { $s = $set; }

        if ($result{'N0'} and not $result{'origin'}) { # skip in case distance is used and no distance is available
            push @{$fpN{$s}{'0'}{'data'}}, $result{'N0'};
            push @{$fpN{$s}{'1'}{'data'}}, $result{'N1'};
            push @{$fpN{$s}{'2'}{'data'}}, $result{'N2'};
        }

        # Show progress
        print "." and $iProgress++ if ($iProtein/$nProtein >= $iProgress/60);
    }
    print "\n";
}
print " done.\n";

## create folder to contain results ##

mkdir "$wkDir/result" unless -d "$wkDir/result";
mkdir "$wkDir/result/statistics" unless -d "$wkDir/result/statistics";

## output raw data for further statistical analysis ##

if ($outRaw) {
    open OUT, ">$wkDir/result/statistics/rawdata.txt";
    print OUT "Query\tSet\tLength\tHits\tSelf\tClose\tDistal\n";
    foreach my $set (sort keys %results) {
        $n = @{$results{$set}};
        for ($i=0; $i<$n; $i++) {
            my %res = %{$results{$set}[$i]};
            next if $res{'origin'};
            next unless $res{'N0'};
            print OUT "$res{'accn'}\t$set\t$res{'length'}\t$res{'n'}\t$res{'N0'}\t$res{'N1'}\t$res{'N2'}\n";
        }
    }
    close OUT;
    print "Raw data are saved in result/statistics/rawdata.txt.\n";
    print "You may conduct further analyses on these data.\n";
    if ($interactive) {
        print "Press Enter to continue, or Ctrl+C to exit:";
        $s = <STDIN>;
    }
}

## graph fingerprints with R ##

if ($graphFp) {
    print "\nGraphing fingerprints with R...";
    $R->startR;
    print "R cannot be started. Make sure it is properly installed in the system.\n" and exit 1 unless $R->is_started();
    if ($plot3D) {
        $R->send("library('rgl')");
        print "\n  You chose to display interactive 3D scatter plots. They will be displayed sequentially. Use mouse to rotate plots. Press Enter in the terminal to move to the next plot.\n";
    }
    foreach my $set (sort keys %fpN) {
        my $fpre = "$wkDir/result/statistics/".("$set." x ($set ne "0")); # prefix of filename
        my $tpost = " of $set" x ($set ne "0"); # postfix of title
        my @gcode = ('Self', 'Close', 'Distal');
        for (0..2) { # send data to R
            @b = @{$fpN{$set}{$_}{'data'}};
            $_ = sprintf("%.3f", $_) for (@b);
            $R->send("x$_<-c(".join (",", @b).")");
            @b = sort{$a<=>$b}@b;
            $R->send("lim$_<-".$b[$#b]); # find proper xlim
            if ($exOutlier) {
                @c = boxplot(@b) if ($exOutlier == 1);
                @c = adjusted_boxplot(@b) if ($exOutlier == 2);
                @c = modified_z(@b) if ($exOutlier == 3);
                if ($b[$#b] > $c[1]) {
                    for ($i=$#b; $i>=0; $i--) {
                        if ($b[$i] <= $c[1]) {
                            $R->send("lim$_<-".$b[$i]);
                            last;
                        }
                    }
                }
            }
        }
        my @xr = ([], [], []); # reference positives
        if (keys %refs) {
            foreach my $key (keys %refs) {
                if ($refs{$key}[0]) {
                    push (@{$xr[$_]}, $refs{$key}[$_+1]) for (0..2);
                }
            }
            for (0..2) {
                @b = @{$xr[$_]};
                $_ = sprintf("%.3f", $_) for (@b);
                $R->send("xr$_<-c(".join (",", @b).")");
            }
        }
        if ($boxPlot) {
            $R->send("pdf('$fpre"."box.pdf',useDingbats=F)");
            $R->send("par(mar=c(4,4,1,1)+0.1,mgp=c(2.5,0.75,0))");
            $R->send("boxplot(x0,x1,x2,names=c('self','close','distal'),xlab='Group',ylab='Weight',main='')");
            $R->send("dev.off()");
        }
        if ($selfLow) { @b = (0, 1, 2); }
        else { @b = (1, 2); }
        if ($histogram) {
            for (@b) {
                $R->send("pdf('$fpre"."hist.".lc($gcode[$_]).".pdf',useDingbats=F)");
                $R->send("par(mar=c(4,4,1,1)+0.1,mgp=c(2.5,0.75,0))");
                $R->send("hist(x$_, breaks=$nBin,freq=F,col='lightgrey',xlab='$gcode[$_] weight',ylab='Probability density',main=''".",xlim=range(0:lim$_)" x ($exOutlier and $exOutlier > 0).")");
                $R->send("dev.off()");
            }
        }
        if ($densityPlot) {
            for (@b) {
                $R->send("pdf('$fpre"."density.".lc($gcode[$_]).".pdf',useDingbats=F)");
                if ($plotRef) { $R->send("par(mar=c(4,4,1,3)+0.1,mgp=c(2.5,0.75,0))"); }
                else { $R->send("par(mar=c(4,4,1,1)+0.1,mgp=c(2.5,0.75,0))"); }
                $R->send("d<-density(x$_".(",bw=bw.nrd0(x$_)*$bwF" x ($bwF and $bwF != 1)).")");
                if ($exOutlier) { $R->send("lim<-range(min(d\$x):lim$_)"); }
                else { $R->send("lim<-range(d\$x[is.finite(d\$x)])"); }
                $R->send("plot(d,lwd=2,xlim=lim,xlab='Weight',ylab='Probability density',main='')");
                if (scalar @{$xr[0]}) {
                    $R->send("par(new=TRUE)");
                    $R->send("plot(density(xr$_".(",bw=bw.nrd0(ab)*$bwF" x ($bwF and $bwF != 1))."),xlim=lim,xaxt='n',yaxt='n',xlab='',ylab='',main='',col=2)");
                    $R->send("axis(4,col=2,col.ticks=2)");
                    $R->send("mtext('Probability density of true positives',side=4,line=-2,col=2)");
                    $R->send("rug(xr$_,ticksize=0.04,lwd=1,col=rgb(1,0,0,0.25))");
                }
                $R->send("dev.off()");
            }
        }
        if ($scatterPlot) {
            if ($selfLow) {
                for (@b) {
                    ($i, $j) = ($_, $_+1);
                    $j = 0 if ($j == 3);
                    $R->send("pdf('$fpre"."scatter.".lc($gcode[$i])."-".lc($gcode[$j]).".pdf',useDingbats=F)");
                    $R->send("par(mar=c(4,4,1,1)+0.1,mgp=c(2.5,0.75,0))");
                    $R->send("plot(x$i,x$j,pch=16,col=rgb(0,0,0,0.25),xlim=range(0:lim$i),ylim=range(0:lim$j),xlab='$gcode[$i] weight',ylab='$gcode[$j] weight')");
                    if (scalar @{$xr[0]}) {
                        $R->send("points(xr$i,xr$j,pch=16,col=rgb(1,0,0,0.5))");
                        $R->send("legend('topright','true positive',pch=16,col='red')");
                    }
                    $R->send("dev.off()");
                }
            } else {
                $R->send("pdf(\"$fpre"."scatter.pdf\",useDingbats=F)");
                $R->send("par(mar=c(4,4,1,1)+0.1,mgp=c(2.5,0.75,0))");
                $R->send("plot(x1,x2,pch=16,col=rgb(0,0,0,0.25),xlim=range(0:lim1),ylim=range(0:lim2),xlab='Close weight',ylab='Distal weight',main='')");
                if (scalar @{$xr[0]}) {
                    $R->send("points(xr1,xr2,pch=16,col=rgb(1,0,0,0.5))");
                    $R->send("legend('topright','true positive',pch=16,col='red')");
                }
                $R->send("dev.off()");
            }
        }
        if ($plot3D) { # This function is not functioning properly
            # if ($^O=~/Win/) { $R->send("windows()"); } elsif ($^O=~/Mac/) { $R->send("quartz()"); } else { $R->send("x11()"); }
            $R->send("plot3d(x0,x1,x2,xlab='Self weight',ylab='Close weight',zlab='Distal weight')");
            print "Displaying 3D plot$tpost. Press Enter to move on.";
            $s = <STDIN>;
        }
    }
    $R->stopR();
    print " done.\n";
    print "Graphs are saved in result/statistics/.\n";
    if ($interactive) {
        print "You may take a look at the graphs before proceeding.\n";
        print "Press Enter to continue, or Ctrl+C to exit:";
        $s = <STDIN>;
    }
}


## Compute the statistics for the whole genome(s), i.e., phyletic pattern, or "fingerprint" ##

print "\nComputing statistics...";
if (($howCO == 4 and ($toolKDE or $toolExtrema)) or ($howCO == 5) or $dipTest) {
    $R->startR;
    print "R cannot be started. Make sure it is properly installed in the system.\n" and exit 1 unless $R->is_started();
    $R->send("library(pastecs)") if ($howCO == 4 and $toolExtrema);
    $R->send("library(diptest)") if $dipTest;
}

foreach my $set (keys %fpN) {
    if ($set eq "0") { print "\n  All protein sets:\n"; } else { print "\n  Protein set $set:\n"; }
    foreach my $key ('0','1','2') {
        next unless @{$fpN{$set}{$key}{'data'}};

        if ($key eq "0") { print "    Self group:\n"; } elsif ($key eq "1") { print "    Close group:\n"; } elsif ($key eq "2") { print "    Distal group:\n"; }

        @a = sort {$a<=>$b} @{$fpN{$set}{$key}{'data'}}; # sort low to high

        my $global_cutoff;
        my $computed_cutoff;
        my $use_global = 0;
        my $half_way = median(@a); ########## if computed cutoff > median, use global cutoff

        # compute basic statistical parameters

        $fpN{$set}{$key}{'n'} = @a;
        $s = 0; $s += $_ for @a;
        $fpN{$set}{$key}{'mean'} = $s/$fpN{$set}{$key}{'n'};
        $s = 0; $s += ($fpN{$set}{$key}{'mean'} - $_)**2 for @a;
        $fpN{$set}{$key}{'stdev'} = sqrt($s/($fpN{$set}{$key}{'n'}-1));
        $fpN{$set}{$key}{'min'} = $a[0];
        $fpN{$set}{$key}{'max'} = $a[$#a];
        $fpN{$set}{$key}{'median'} = median(@a);
        $fpN{$set}{$key}{'mad'} = mad(@a);
        ($fpN{$set}{$key}{'q1'}, $fpN{$set}{$key}{'q3'}) = quantiles(@a);

        if ($key eq '0' and not $selfLow) {
            $fpN{$set}{$key}{'cutoff'} = 0;
            print "      Skipped.\n";
            next;
        }

        # determine cutoff using global cutoff

        $i = $fpN{$set}{$key}{'n'}*$globalCO;
        if (int($i) == $i) {
            $global_cutoff = ($a[$i-1]+$a[$i])/2;
        } else {
            if ($i-int($i) <= int($i)-$i+1) {
                $global_cutoff = $a[int($i)];
            } else {
                $global_cutoff = $a[int($i)-1];
            }
        }
        print "      Global cutoff ($globalCO) = $global_cutoff.\n";

        # override individual cutoff with global cutoff

        $use_global = 1 if ((($key eq '0') and ($selfCO eq 'G')) or (($key eq '1') and ($closeCO eq 'G')) or (($key eq '2') and ($distalCO eq 'G')));

        # exclude outliers

        if ($exOutlier) {
            @c = boxplot(@a) if ($exOutlier == 1);
            @c = adjusted_boxplot(@a) if ($exOutlier == 2);
            @c = modified_z(@a) if ($exOutlier == 3);
            if ($a[$#a] > $c[1]) {
                for ($i=$#a; $i>=0; $i--) {
                    if ($a[$i] <= $c[1]) {
                        @a = @a[0..$i];
                        last;
                    }
                }
            }
        }

        # perform Hartigan's dip test to assess non-unimodality

        if ($dipTest) {
            print "      Performing Hartigan's dip test...";
            $R->send("x<-c(".join (",", @a).")");
            $R->send("dip.test(x)");
            $s = $R->read;
            open OUT, ">$key.out";
            print OUT join(' ', @a);
            close OUT;

            if ($s =~ /D = (\S+), p-value [<=] (\S+)\n/) {
                print " done.\n";
                print "      D = $1, p-value = $2\n";
                if ($dipSig) {
                    if ($2 >= $dipSig) {
                        print "      The weight distribution is NOT significantly non-unimodal.\n";
                        if ($howCO >= 3) {
                            if ($interactive) {
                                print "      Proceed with statistical analysis anyway (yes) or use global cutoff ($globalCO) instead (NO)? ";
                                while (<STDIN>) {
                                    chomp;
                                    unless ($_) { $use_global = 1; last; }
                                    if (/^y$/i or /^yes$/i) { last; }
                                    if (/^n$/i or /^no$/i) { $use_global = 1; last; }
                                }
                            } else { $use_global = 1; }
                        }
                    } else { print "      The weight distribution is significantly non-unimodal.\n"; }
                } else { print "      The weight distribution is ". ("NOT " x ($2 >= 0.05)). "significantly non-unimodal.\n"; }
            } else { print " failed.\n"; }
        }

        # determine cutoff using histogram

        print "      Generating histogram..." if ($howCO == 3 and not $use_global);

        my @freqs = (0)x$nBin;
        my $interval = ($a[$#a]-$a[0])/$nBin;
        my $cid = 0;
        for ($j=1; $j<$nBin; $j++) {
            my $high_bound = $interval*$j;
            for ($i=$cid; $i<=$#a; $i++) {
                if ($a[$i] < $high_bound) {
                    $freqs[$j-1] ++;
                } else {
                    $cid = $i;
                    last;
                }
            }
        }
        $freqs[$nBin-1] = $#a-$cid+1;
        my $local_min = 0; # index of the lowest bar from left
        for ($i=1; $i<$nBin-1; $i++) {
            if ($freqs[$i]<$freqs[$i-1] and $freqs[$i]<=$freqs[$i+1]) {
                $computed_cutoff = $interval*$i;
                $local_min = $i;
                last;
            }
        }

        print " done.\n" if ($howCO == 3 and not $use_global);

        # draw histogram

        if ($plotHist) {
            print "  Histogram:";
            @c = sort {$b<=>$a} @freqs;
            $s = 50/$c[0];
            my @widths = (0)x$nBin;
            my @labels = (0)x$nBin;
            for ($i=0; $i<$nBin; $i++) {
                $widths[$i] = int($freqs[$i]*$s);
                $labels[$i] = sprintf("%.2f", $interval*$i)."-".sprintf("%.2f", $interval*($i+1));
            }
            @c = sort {length($b)<=>length($a)} @labels;
            $s = length($c[0]);
            print "\n";
            for ($i=0; $i<$nBin; $i++) {
                print " "x($s-length($labels[$i]))."$labels[$i] ".("*"x$widths[$i])."$freqs[$i]".(" (local minimum)"x($local_min and $local_min==$i))."\n";
            }
        }

        # determine cutoff using kernel density estimation (KDE)

        if ($howCO == 4 and not $use_global) {

            # perform kernel density estimation (KDE)

            print "      Performing kernel density estimation...";
            my (@dx, @dy); # x- and y-coornidates of density function
            my $KDEstats = "";

            # perform KDE using basic R command "density"

            if ($toolKDE) {
                $R->send("x<-c(".join (",", @a).")");
                if ($bwF and $bwF != 1) {
                    $R->send("bwx<-bw.nrd0(x)*$bwF");
                    $R->send("d<-density(x,bw=bwx)");
                } else {
                    $R->send("d<-density(x)");
                }
                @dx = @{$R->get('d$x')};
                @dy = @{$R->get('d$y')};
            }

            # perform KDE using self-written Perl code

            else {
                my $n = scalar @a;
                my $mean; $mean += $_ for @a; $mean = $mean/$n;
                my $stdev; $stdev += ($mean-$_)**2 for @a; $stdev = sqrt($stdev/($n-1));
                my @Q = quantiles(@a); my $iqr = $Q[1]-$Q[0];
                my $bw;
                if ($stdev == 0 and $iqr == 0) { $bw = 1; }
                elsif ($stdev == 0) { $bw = $iqr/1.34; }
                elsif ($iqr == 0) { $bw = $stdev; }
                elsif ($stdev <= $iqr/1.34) { $bw = $stdev; }
                else { $bw = $iqr/1.34; }
                $bw = 0.9*$bw*$n**(-1/5); # select bandwidth by Silverman's ¡°rule of thumb¡± (1986)
                $bw = $bw*$bwF if ($bwF and $bwF != 1);
                my ($min, $max) = ($a[0]-3*$bw, $a[$#a]+3*$bw); # cut = 3
                $KDEstats = "      N = $n, bandwidth = ".sprintf("%.3f", $bw).".\n";
                for (my $x=$min; $x<=$max; $x+=($max-$min)/511) { # 512 points
                    my $e = 0; $e += exp(-(($x-$_)/$bw)**2/2)/sqrt(2*3.1415926536) for @a; # Gaussian kernel
                    push @dx, $x;
                    push @dy, 1/$n/$bw*$e;
                }
            }

            print " done.\n";
            print $KDEstats;

            # plot density function on screen

            if ($plotKDE) {
                print "  Density function:\n";
                my $k_x = int(@dx/100);
                $k_x = int(@dx/$plotKDE) if ($plotKDE > 1);
                my $k_y = 0;
                foreach (@dy) {
                    $k_y = $_ if ($_ > $k_y);
                }
                $k_y = 64/$k_y;
                for ($i=0; $i<=$#dy; $i+=$k_x) {
                    print " "x(7-length(sprintf("%.2f", $dx[$i]))).sprintf("%.2f", $dx[$i]);
                    print " "x (int($dy[$i]*$k_y)+1)."*\n";
                }
            }

            my $peak1st;
            my $failed;
            my ($peak_i, $peak_x, $peak_y);
            my ($pit_i, $pit_x, $pit_y);

            # identify local extrema using R package "pastecs"

            if ($toolExtrema) {
                if ($toolKDE < 2) {
                    $R->send("d<-data.frame(x=c(".join (",", @dx)."),y=c(".join (",", @dy)."))");
                }
                $R->send("tp<-turnpoints(ts(d\$y))");
                if ($R->get('tp$firstispeak') eq "TRUE") {
                    if ($R->get('tp$nturns') == 1) { $failed = 1; }
                    else { $peak1st = 1; }
                } else {
                    $peak1st = 0;
                }
                unless ($failed) {
                    my @tpx = @{$R->get("d\$x[tp\$tppos]")};
                    my @tpy = @{$R->get("d\$y[tp\$tppos]")};
                    if ($peak1st) {
                        ($peak_x, $peak_y) = ($tpx[0], $tpy[0]);
                        ($pit_x, $pit_y) = ($tpx[1], $tpy[1]);
                    } else {
                        ($peak_x, $peak_y) = ($dx[0], $dy[0]);
                        ($pit_x, $pit_y) = ($tpx[0], $tpy[0]);
                    }
                }
            }

            # identify local extrema using self-written Perl code

            else {
                for ($i=0; $i<$#dy; $i++) {
                    if ($dy[$i] < $dy[$i+1]) { $peak1st = 1; last; }
                    elsif ($dy[$i] > $dy[$i+1]) { $peak1st = 0; last; }
                    else { next; }
                }
                if ($peak1st) {
                    for ($i=1; $i<$#dy; $i++) {
                        if (($dy[$i-1] <= $dy[$i]) and ($dy[$i] > $dy[$i+1])) {
                            ($peak_i, $peak_x, $peak_y) = ($i, $dx[$i], $dy[$i]);
                        }
                        if (($dy[$i-1] > $dy[$i]) and ($dy[$i] <= $dy[$i+1])) {
                            ($pit_i, $pit_x, $pit_y) = ($i, $dx[$i], $dy[$i]);
                        }
                        last if ($peak_i and $pit_i);
                    }
                    $failed = 1 unless $pit_i;
                } else {
                    ($peak_i, $peak_x, $peak_y) = (0, $dx[0], $dy[0]);
                    for (my $i=1; $i<$#dx; $i++) {
                        if (($dy[$i-1] > $dy[$i]) and ($dy[$i] <= $dy[$i+1])) {
                            ($pit_i, $pit_x, $pit_y) = ($i, $dx[$i], $dy[$i]);
                        }
                        last if $pit_i;
                    }
                }
            }

            if ($failed) {
                print "      The weight population cannot be clustered using kernel density estimation. This is typically caused by an even distribution of all weights.\n";
                $use_global = 1;
            } else {

                # locate cutoff point

                if ($modKCO == 0) { # pit
                    $computed_cutoff = $pit_x;
                } elsif ($modKCO == 1) { # horizontal midpoint
                    $computed_cutoff = ($pit_x+$peak_x)/2;
                } elsif ($modKCO == 2) { # horizontal quantile
                    $computed_cutoff = $pit_x-($pit_x-$peak_x)*$qKCO;
                } else { # vertical quantile
                    my $vCO = 0;
                    $vCO = $pit_y+($peak_y-$pit_y)*$qKCO;
                    for (my $i=$peak_i; $i<$pit_i; $i++) {
                        if (($dy[$i] >= $vCO) and ($vCO > $dy[$i+1])) {
                            $computed_cutoff = $dx[$i]+($dx[$i+1]-$dx[$i])*($dy[$i]-$vCO)/($dy[$i]-$dy[$i+1]);
                            last;
                        }
                    }
                }
            }
        }

        # determine cutoff using hierarchical clustering

        if ($howCO == 5 and not $use_global) {
            print "      Performing hierarchical clustering...";
            $R->send("x<-c(".join (",", @a).")");
            $R->send("d<-dist(x,method=\"euclidean\")");
            $R->send("fit<-hclust(d,method=\"ward\")");
            print " done.\n";
            my $nCluster = 1;
            while ($nCluster++) {
                $R->send("c<-cutree(fit,k=$nCluster)");
                my $clusters = $R->get('c'); # tip: return vector from R
                @c = (()) x $nCluster; # data of clusters
                for ($i=0; $i<=@{$clusters}-1; $i++) {
                    foreach (1..$nCluster) {
                        if (@{$clusters}[$i] == $_) {
                            push @{$c[$_-1]}, $a[$i];
                            last;
                        }
                    }
                }
                for ($i=1; $i<=$#c; $i++) {
                    @{$c[$i]} = sort {$a <=> $b} @{$c[$i]};
                }
                @b = sort {$a->[0] <=> $b->[0]} @c;
                print "      With $nCluster clusters, cutoff is $b[1][0]. Accept? (YES/no) ";
                my $user_okay = 1;
                while (<STDIN>) {
                    chomp; last unless $_;
                    last if (/^y$/i or /^yes$/i);
                    if (/^n$/i or /^no$/i) { $user_okay = 0; last; }
                }
                last if $user_okay;
            }
            $computed_cutoff = $b[1][0];
        }

        # report cutoff to user

        if ($howCO == 0 or ($howCO >= 3 and $use_global)) { # user-defined global cutoff
            $fpN{$set}{$key}{'cutoff'} = $global_cutoff;
            print "      Cutoff is $global_cutoff (determined by global cutoff $globalCO)\n";
        }
        if ($howCO == 1 and $unite) { # user-defined individual cutoff
            if ($key eq '0') { $s = $selfCO; }
            elsif ($key eq '1') { $s = $closeCO; }
            elsif ($key eq '2') { $s = $distalCO; }
            if ($s) {
                $fpN{$set}{$key}{'cutoff'} = $s;
                print "      Cutoff is $s (user-defined).\n";
            } else {
                $fpN{$set}{$key}{'cutoff'} = $global_cutoff;
                print "      User-defined cutoff is not available. Use global cutoff $global_cutoff instead.\n";
            }
        }
        if ($howCO >= 3 and not $use_global) { # computed cutoff
            $s = 0; if ($howCO == 3) { $s = "histogram"; }
            elsif ($howCO == 4) { $s = "kernel density estimation"; }
            elsif ($howCO == 5) { $s = "clustering analysis"; }
            if ($computed_cutoff) {
                if ($computed_cutoff <= $half_way) {
                    $fpN{$set}{$key}{'cutoff'} = $computed_cutoff;
                    print "      Cutoff is ".sprintf("%.3f", $computed_cutoff)." (determined by $s).\n";
                } else {
                    $fpN{$set}{$key}{'cutoff'} = $global_cutoff;
                    print "      ". ucfirst($s) ." identified a cutoff ".sprintf("%.3f", $computed_cutoff)." which is too large. Use global cutoff $global_cutoff instead.\n";
                }
            } else {
                $fpN{$set}{$key}{'cutoff'} = $global_cutoff;
                print "      ". ucfirst($s) ." failed to identify a cutoff. Use global cutoff $global_cutoff instead.\n";
            }
        }
        # $fpN{$set}{$key}{'cutoff'} = 0.00001 unless $fpN{$set}{$key}{'cutoff'};

        # ask user to enter cutoff

        my $user_enter;
        if ($interactive and $fpN{$set}{$key}{'cutoff'}) {
            print "      Accept? (YES/no) ";
            while (<STDIN>) {
                chomp;
                last unless $_;
                last if (/^y$/i or /^yes$/i);
                if (/^n$/i or /^no$/i) { $user_enter = 1; last; }
            }
        }
        if ($user_enter or $howCO == 2) {
            print "      Enter user-specified cutoff: ";
            while (<STDIN>) {
                chomp;
                last if /^\d+\.?\d*$/;
                print "      Invalid cutoff value. Re-enter: ";
            }
            $fpN{$set}{$key}{'cutoff'} = $_;
        }
    }
}
if (($howCO == 4 and ($toolKDE or $toolExtrema)) or ($howCO == 5) or $dipTest) {
    $R->stopR;
}
print " done.\n";


## output fingerprint ##

if ($outFp) {
    open OUT, ">$wkDir/result/statistics/fingerprint.txt";
    print OUT "#NEXUS\nBEGIN STATISTICS;\n";
    print OUT "\tGroup\tNumber\tMean\tSD\tMin\tMax\tMedian\tMAD\tQ1\tQ3\tCutoff\n";
    foreach my $set (sort keys %fpN) {
        foreach ('0','1','2') {
            %h = %{$fpN{$set}{$_}};
            if ($_ eq '0') { $s = "self"; } elsif ($_ eq '1') { $s = "close"; } elsif ($_ eq '2') { $s = "distal"; }
            if ($set eq '0') { print OUT "all"; } else { print OUT $set; }
            print OUT "\t$s\t$h{'n'}\t".sprintf("%.2f", $h{'mean'})."\t".sprintf("%.2f", $h{'stdev'})."\t".sprintf("%.2f", $h{'min'})."\t".sprintf("%.2f", $h{'max'})."\t".sprintf("%.2f", $h{'median'})."\t".sprintf("%.2f", $h{'mad'})."\t".sprintf("%.2f", $h{'q1'})."\t".sprintf("%.2f", $h{'q3'})."\t".sprintf("%.2f", $h{'cutoff'})."\n";
        }
    }
    print OUT "END;\n";
    close OUT;
    print "Result is saved in result/statistics/fingerprint.txt.\n";
}

if ($interactive) {
    print "Press Enter to proceed with prediction, or Ctrl+C to exit:";
    $s = <STDIN>;
}


## conduct bidirectional best hit (BBH) search ##

if ($BBH == 2) {
    print "Conducting bidirectional best hit (BBH) search...\n";
    unless (-e "$wkDir/result/bbh_input.txt") {
        open OUT, ">$wkDir/result/bbh_input.txt";
        foreach my $set (keys %results) {
            for ($i=0; $i<@{$results{$set}}; $i++) {
                my %res = %{$results{$set}[$i]};
                if (exists $res{'BBH'} and $res{'BBH'}) {
                    # set id query_accn query_taxid subject_accn subject_taxid
                    print OUT $set."\t".$i."\t".$res{'accn'}."\t".$selfinfo{$set}{'taxid'}."\t".$res{'hit2accn'}."\t".$res{'hit2'}."\n";
                }
            }
        }
    }
    $s = $0; $s =~ s/analyzer\.pl$/bbh.pl/;
    system "$^X $s $wkDir";
    unlink "$wkDir/result/bbh_input.txt";
    open IN, "<$wkDir/result/bbh.txt";
    while (<IN>) {
        s/\s+$//;
        @a = split (/\t/);
        $results{$a[0]}[$a[1]]{'BBH'} = "" if ($a[$#a] ne '1');
    }
    close IN;
    print " done.\n";
}


################################################
## Analyze the whole master record and output ##
################################################

mkdir "$wkDir/result/detail" unless -d "$wkDir/result/detail";

print "Predicting...";
foreach my $set (keys %results) {
    $n = @{$results{$set}};
    for ($i=0; $i<$n; $i++) {
        my %res = %{$results{$set}[$i]};
        next if ($minHits && ($res{'n'} < $minHits));
        next if ($minSize && ($res{'length'} < $minSize));

        my $fp = $set;     $fp = '0' if $unite;

        # principle: fewer hits suggests income or loss
        # check if self hits are low (indicating the gene is not prevalent within the group)
        if (($res{'N0'} < $fpN{$fp}{'0'}{'cutoff'}) or ($res{'N0'} == 0)) {
            $res{'in'} = 1;
        }

        ## predict incoming HGT events ##
        # hits from close sister groups are low, indicating there's no vertical ancestor
        if (($res{'N1'} < $fpN{$fp}{'1'}{'cutoff'}) or ($res{'N1'} == 0)) {
            # overall non-self hits are normal, indicating it's not an origination event
            if ($res{'N2'} and $res{'N2'} >= $fpN{$fp}{'2'}{'cutoff'}) {
                $res{'income'} = 1;
                $res{'income'} = "" if ($selfLow and (($res{'N0'} > 0) and ($res{'N0'} >= $fpN{$fp}{'0'}{'cutoff'})));
            }
        }

        ## predict gene loss events ##
        # hits from close sister groups are normal,
        elsif ($res{'in'}) {
            # self hits are low, suggesting the gene is absent in some lineages
            $res{'loss'} = 1;
        }

        if (0) { # not available in this version
            ## predict gene origination event ##
            # detection won't work for saturated search result.
            if ($maxHits and ($res{'N0'}+$res{'N1'}+$res{'N2'}) < $maxHits) {
                # if every hit is self, then must be origination.
                unless ($res{'N1'}+$res{'N2'}) {
                    $res{'origin'} = 1;
                } else {
                    # if non-self hits are significantly weak.
                    $s = 0;
                    $s = $res{'S1'}[0] if $res{'N1'};
                    $s = $res{'S2'}[0] if ($res{'N2'} and $res{'S2'}[0] > $s);
                }
            }
        }

        # predicting gene outcoming events
        # self hits must be normal, excluding any paralogs
        unless ($res{'in'}) {
            if ($res{'N1'}+$res{'N2'} and exists($res{'sGenus'}{'q1'})) {
                if ($res{'nGenus'}{'max'} > $res{'sGenus'}{'q1'}) {        # the best hit from other genera fall within the range of self genus hits
                    @a = @{$res{'nGenus'}{'scores'}};
                    $s = 0;
                    foreach (@a) {
                        if ($_ < $res{'sGenus'}{'min'}) { last;}
                        else { $s = $_;}
                    }
                    if ($s > $res{'sGenus'}{'q1'}) {
                        $res{'genus_outcome'} = $res{'nGenusHit'};
                    }
                }
            }
        }
        $results{$set}[$i] = {%res};
    }

    ###### Output report ######

    open (OUT, ">$wkDir/result/detail/$set.txt");
    @a = ('Query','Length','Product','Hits','Self','Close','Distal','HGT');
    unless ($BBH) {
        push (@a, "Loss") if $loss;
        push (@a, "POE") if $POE;
    }
    push @a, "Match";
    print OUT "HGTector result of $set\n".join("\t",@a)."\n";
    for ($i=0; $i<$n; $i++) {
        %h = %{$results{$set}[$i]};
        print OUT $h{'query'}."\t".$h{'length'}."\t".$h{'product'}."\t".$h{'n'}."\t";
        print OUT "\n" and next if ($minHits and ($h{'n'} < $minHits)) or ($minSize and ($h{'length'} < $minSize)); #####??????#######
        print OUT sprintf("%.2f", $h{'N0'})."\t".sprintf("%.2f", $h{'N1'})."\t".sprintf("%.2f", $h{'N2'})."\t";
        unless ($BBH) {
            print OUT $h{'income'};
            print OUT "\t".$h{'loss'} if $loss;
            print OUT "\t".$h{'origin'} if $POE;
        } else {
            print OUT $h{'BBH'} if exists $h{'BBH'} and $h{'BBH'};
        }
        print OUT "\t";
        if (exists $h{'hit2'} and exists $taxadb{$h{'hit2'}}) {
            print OUT $h{'hit2'}." (".$taxadb{$h{'hit2'}}{'name'}.")";
        }
        print OUT "\n";
    }
    close OUT;
}
print " done.\n";
print "Prediction results are saved in result/detail/.\n";
exit 0;



##################
## sub routines ##
##################

sub median (@) {
    my $mid = int(@_/2);
    return $_[0] unless $mid;
    if (@_ % 2) {
        return $_[$mid];
    } else {
        return ($_[$mid-1] + $_[$mid])/2;
    }
}

sub mad (@) {
    my @absdev;
    my $median = median(@_);
    foreach my $x(@_) {
        push @absdev, abs($x - $median);
    }
    return median(sort{$a<=>$b}@absdev);
}

sub quantiles (@) {
    my $Q1, my $Q3;
    my $mid = int(@_/2);
    return ($_[0],$_[0]) unless $mid;
    if (@_ % 2) {
        $Q1 = median(@_[0..$mid]);
        $Q3 = median(@_[$mid..$#_]);
    } else {
        $Q1 = median(@_[0..($mid-1)]);
        $Q3 = median(@_[$mid..$#_]);
    }
    return ($Q1,$Q3);
}


# compute Z-score (Z = (xi-x^)/s)

sub z_scores(@) {
    my $mean = 0;
    $mean += $_ for @_;
    $mean = $mean / @_;
    my $stdev = 0;
    $stdev += ($mean-$_)**2 for @_;
    $stdev = sqrt($stdev/(@_-1));
    my @z = ();
    push (@z, ($_-$mean)/$stdev) for @_;
    return @z;
}


# Z-score test for outliers (|Z| > 3)

sub z_test(@) {
    my @data = sort{$a<=>$b}@_;
    my @z = z_scores(@data);
    my $lower_fence = $data[0];
    my $upper_fence = $data[$#data];
    for (my $i=0; $i<=$#data; $i++) {
        if (abs($z[$i]) <= 3) {
            $lower_fence = $data[$i];
            last;
        }
    }
    for (my $i=$#data; $i>=0; $i--) {
        if (abs($z[$i]) <= 3) {
            $upper_fence = $data[$i];
            last;
        }
    }
    return ($lower_fence, $upper_fence);
}

# modified Z-score test for outliers (|modified_Z| > 3.5) (Iglewicz and Hoaglin, 1993)

sub modified_z(@) {
    my @data = sort{$a<=>$b}@_;
    my $lower_fence = $data[0];
    my $upper_fence = $data[$#data];
    my $median = median(@data);
    my $mad = mad(@data);
    return ($data[0],$data[$#data]) unless $mad;
    for (my $i=0; $i<=$#data; $i++) {
        if (abs(0.6745*($data[$i]-$median)/$mad) <= 3.5) {
            $lower_fence = $data[$i];
            last;
        }
    }
    for (my $i=$#data; $i>=0; $i--) {
        if (abs(0.6745*($data[$i]-$median)/$mad) <= 3.5) {
            $upper_fence = $data[$i];
            last;
        }
    }
    return ($lower_fence, $upper_fence);
}


# boxplot test for outliers

sub boxplot(@) {
    my $lower_fence, my $upper_fence;
    my @data = sort{$a<=>$b}@_;
    my @Q = quantiles(@data);
    my $iqr = $Q[1]-$Q[0];
    my $f = 3*exp(10/@data);
    $lower_fence = $Q[0]-$f*$iqr;
    $upper_fence = $Q[1]+$f*$iqr;
    return ($lower_fence, $upper_fence);
}


# adjusted boxplot test for outliers

sub adjusted_boxplot(@) {
    my $lower_fence, my $upper_fence;
    my @data = sort{$a<=>$b}@_;
    my $median = median(@data);
    my @lower; my @upper;
    foreach (@data) {
        push (@lower, $_) if ($_ <= $median);
        push (@upper, $_) if ($_ >= $median);
    }
    my @kernel;
    foreach my $i (@lower) {
        foreach my $j (@upper) {
            next if ($i == $j);
            push @kernel, (($j-$median)-($median-$i))/($j-$i);
        }
    }
    my $mc = median(sort{$a<=>$b}@kernel);
    my @Q = quantiles(@data);
    my $iqr = $Q[1]-$Q[0];
    my $f = 1.5; # *exp(10/@data);
    if ($mc >= 0) {
        $lower_fence = $Q[0]-$f*exp(-3.5*$mc)*$iqr;
        $upper_fence = $Q[1]+$f*exp(4*$mc)*$iqr;
    } else {
        $lower_fence = $Q[0]-$f*exp(-4*$mc)*$iqr;
        $upper_fence = $Q[1]+$f*exp(3.5*$mc)*$iqr;
    }
    return ($lower_fence, $upper_fence);
}

sub recurse_deOutlier(@) { # assume data are high to low
    my @data = @_;
    my @fences = adjusted_boxplot @data;
    for ($i=0; $i<=$#data; $i++) {
        if ($data[$i] < $fences[0]) {
            $i--;
            last;
        }
    }
    if ($i < $#data) {
        @data = @data[0..$i];
        @data = recurse_deOutlier @data;
    }
    return @data;
}

sub recurse_Z(@) {
    my @data = @_;
    if (modified_z(@data) > 3.5) {
        pop @data;
        @data = recurse_Z(@data);
    }
    return @data;
}
