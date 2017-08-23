#!/usr/bin/perl

use warnings;
use strict;
$| = 1;

print "
-> Orthologer: Identify orthologous groups based on homology search results. <-
";

print "

Usage:
  perl orthologer.pl <working directory>

Output:
  taxonomy/orthology.db

" and exit unless @ARGV;


## Global variables ##

my $i; my $j; my $n; my $s; my @a; my @b; my @c; my %h;


## Program parameters ##

my $wkDir = $ARGV[0];               # working directory
my @sets = ();                      # list of protein sets

my %taxa = ();                      # set -> taxid
my %proteins = ();                  # all input proteins (accn -> id)

my @cogs = ();                      # all COGs. array of hashes, including:
                                    # accns -> accn1/accn2/accn3...
                                    # names -> name1    name2    name3...
                                    # nProtein -> number of proteins
                                    # nSet -> number of genomes
                                    # name -> best name
my %names = ();                     # id -> names
my %accns = ();                     # id -> accns
my $cogid = 1;                      # current id

my $evalue = 1e-5;                  # maximum E-value threshold

my $pairRule = 3;                   # pairing rule

## Read configuration ##

if (-e "$wkDir/config.txt") {
    open IN, "<$wkDir/config.txt";
    while (<IN>) {
        s/#.*$//; s/\s+$//g; s/^\s+//g; next unless $_;
        @sets = split (/,/, $1) if /^inSets=(.+)$/;
        $evalue = $1 if /^evalue=(.+)$/;
        $pairRule = $1 if /^pairRule=(\d)$/;
    }
}

unless (@sets) {
    opendir (DIR, "$wkDir/search");
    @a = readdir(DIR);
    close DIR;
    foreach (@a) {
        next if /^\./;
        next unless -d "$wkDir/search/$_";
        push @sets, $_;
    }
}

## Read input proteins ##

print "Identifying orthologous groups (OGs)...\n";
print "  Reading input proteins... ";
foreach my $set (@sets) {
    next unless -d "$wkDir/search/$set";
    opendir (DIR, "$wkDir/search/$set");
    my @blasts = grep(/\.txt$/,readdir(DIR));
    close DIR;
    foreach (@blasts) {
        /^(.+)\.txt$/;
        %h = ('set', $set, 'id', "");
        $proteins{$1} = {%h};
    }
}
print " done. ".(scalar keys %proteins)." proteins read.\n";


## Read input taxonomy ##
open IN, "<$wkDir/taxonomy/self.info";
while (<IN>) {
    @a = split(/\t/);
    $taxa{$a[0]} = $a[2];
}
close IN;


##########################################################
## Parsing BLAST results and identify orthologous pairs ##
##########################################################

print "  Parsing BLAST results...\n    ";
foreach my $set (@sets) {
    next unless -d "$wkDir/search/$set";
    opendir (DIR, "$wkDir/search/$set");
    my @blasts = grep(/\.txt$/,readdir(DIR));
    close DIR;
    my $myTaxon = $taxa{$set};
    print "$set ";
    foreach (@blasts) {
        /^(.+)\.txt$/;
        my $id; # array id of OGs
        my $myAccn = $1;
        next unless exists $proteins{$myAccn};
        open IN, "<$wkDir/search/$set/$_" or next;
        my $product;
        my $reading = 0;
        @b = (); # all accns occurred in this report
        %h = (); # store used taxids, to rule out duplicates.
        while (<IN>) {
            s/\s+$//; next unless $_;
            last if $reading and /^;/;
            if (/^\tProduct=(.+);$/) { $product = $1; $product =~ s/\s+$//; }
            if (/BEGIN ORGANISM;/) { $reading = 1; next; }
            next unless $reading;
            @a = split (/\t/);
            next if $a[4] > $evalue;
            if ($pairRule == 2 or $pairRule == 4 or $pairRule == 5) { # must be best hit from another organism
                next if $a[2] eq $myTaxon;
                next if exists $h{$a[2]};
                $h{$a[2]} = 1;
            }
            @a = split (/\//, $a[0]);
            foreach (@a) {
                next if $_ eq $myAccn;
                if (exists $proteins{$_}) {
                    next if $proteins{$_}{'id'}; # if a protein is already assigned an OG ID, then skip
                    push (@b, $_);
                }
            }
        }
        close IN;

        if ($proteins{$myAccn}{'id'}) { # existed OG
            $id = $proteins{$myAccn}{'id'};
        } else { # new OG
            %h = ('accns', $myAccn, 'names', $product);
            push @cogs, {%h};
            $id = $#cogs;
            $proteins{$myAccn}{'id'} = $id;
        }

        # check reversal Blast results

        foreach my $accn (@b) {
            next if $proteins{$accn}{'id'};
            $product = "";
            $reading = 0;
            my $found = 0; # the source gene is found in the reverse Blast report
            open IN, "<$wkDir/search/$proteins{$accn}{'set'}/$accn.txt" or next;
            while (<IN>) {
                s/\s+$//; next unless $_;
                last if $reading and /^;/;
                if (/^\tProduct=(.+);$/) {
                    $product = $1;
                    $product =~ s/\s+$//;
                    last if ($pairRule < 3); # single direction
                }
                if (/BEGIN ORGANISM;/) { $reading = 1; next; }
                next unless $reading;
                @a = split (/\t/);
                next if $a[4] > $evalue;
                my $taxon = $a[2];
                @a = split (/\//, $a[0]);
                foreach (@a) {
                    if ($myAccn eq $_) {
                        $found = 1;
                        last;
                    }
                }
                last if $found;
                last if (($pairRule == 4) and ($myTaxon eq $taxon) and !$found);
            }
            close IN;
            if ($pairRule >= 3) { # bidirectional hits
                unless ($found) {
                    # print "not match! $set $myAccn <> $accn.\n";
                    next;
                }
            }
            $proteins{$accn}{'id'} = $id;
            $cogs[$id]{'names'} .= "\t$product" if $product;
            $cogs[$id]{'accns'} .= "/$accn";
        }
    }
}
print "\n  done. ".(scalar @cogs)." OGs identified.\n";


## Find out a best name for each COG ##
# rule: the shortest name without "hypothetical"

print "  Naming OGs...";
for ($i=0; $i<=$#cogs; $i++) {
    @a = split(/\t/, $cogs[$i]{'names'});
    @a = sort {length($a) <=> length($b)}(@a);
    my @noHypo = ();
    foreach (@a) {
        push (@noHypo, $_) unless (/hypothetical/ or /hypotethical/ or /hypothetcial/);
    }
    if (@noHypo) { $cogs[$i]{'name'} = "$noHypo[0]"; }
    else{ $cogs[$i]{'name'} = "$a[0]"; }
}
print " done.\n";


#    my $shortName; my $longName;
#    my @noHypo; # subset of names which does not contain "hypothetical"
#    foreach (@a) {
#        push (@noHypo, $_) unless (/hypothetical/ or /hypotethical/ or /hypothetcial/);
#    }
#    if (@noHypo) {
#        foreach (@noHypo) {
#            if (/^([A-Za-z0-9]+) gene product$/ or /^([A-Za-z0-9]+) protein$/ or /^protein ([A-Za-z0-9]+)$/) {
#                $shortName = $1; last;
#            }
#        }
#        @noHypo = sort {length($a) <=> length($b)}(@noHypo);
#        $longName = $noHypo[$#noHypo];
#        $shortName or $shortName = $noHypo[0];
#    } else {
#        @a = sort {length($a) <=> length($b)}(@a);
#        $longName = $a[$#a];
#        $shortName or $shortName = $a[0];
#    }
#    $names{$key} = $shortName;
#    # $names{$key} = "$shortName\t$longName";


## Sort COGs ##

print "  Sorting OGs...";
for ($i=0; $i<=$#cogs; $i++) {
    @a = split(/\//, $cogs[$i]{'accns'});
    %h = ();
    $h{$proteins{$_}{'set'}} = 1 for @a;
    $cogs[$i]{'nProtein'} = scalar (@a);
    $cogs[$i]{'nSet'} = scalar (keys %h);
}
print " done.\n";

@cogs = sort {$b->{'nProtein'} <=> $a->{'nProtein'}} @cogs;


## Export COGs ##

open OUT, ">$wkDir/taxonomy/orthology.db";
for ($i=0; $i<=$#cogs; $i++) {
    print OUT ($i+1)."|$cogs[$i]{'name'}\t$cogs[$i]{'accns'}\n";
}
close OUT;

print "  OG names and members exported as taxonomy/orthology.db.\n";

exit 0;
