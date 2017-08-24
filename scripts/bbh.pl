#!/usr/bin/perl

use strict;
use warnings;
$| = 1;

print "

Performs reversal BLAST, as part of the bidirectional best hit approach.

Usage:
  perl bbh.pl <working directory>

Output:
  Will be appended to the input files.

" and exit unless @ARGV;

print "Performing reversal BLAST...\n";


## public variables ##

my $i; my $j; my $s; my $t; my @a; my @b; my %h;

## global variables ##

my @queries = ();                   # list of query proteins
my %done = ();                      # queries that is already done.

my $query;                          # current query protein (accession number)
my $taxid;                          # target organism
my $subject;                        # subject protein (accession number)
my $eqText;                         # entrez query text
my $isBBH;                          # whether a hit is BBH

my $iRetry = 0;                     # current number of retries
my $isError = 0;                    # if error

use LWP::UserAgent;
use HTTP::Request::Common qw(POST);
my $ua = LWP::UserAgent->new;

## subroutines ##

sub blast ();
sub retry ();

## program parameters ##

my $wkDir = $ARGV[0];

my $httpBlast = 0;                                    # BLAST mode (0: http, 1: local, 2: remote)
my $maxHits = 0;                                    # maximum number of valid hits to preserve. if 0 then = nHit

my $retries = 10;                                    # maximum number of retries
my $nHit = 100;                                        # number of hits to return
my $evalue = 0.01;                                    # maximum E-value cutoff

my $blastServer = "http://blast.ncbi.nlm.nih.gov/Blast.cgi";
my $protdb = "nr";

my $blastdbcmd = "blastdbcmd";
my $blastp = "blastp";
my $threads = 1;                                    # Multiple threads for local BLAST program

## read configurations ##

if (-e "$wkDir/config.txt") {
    open IN, "<$wkDir/config.txt";
    while (<IN>) {
        s/#.*$//; s/\s+$//g; s/^\s+//g; next unless $_;
        $httpBlast = $1 if /^httpBlast=(\d)$/;
        $nHit = $1 if /^nHit=(\d+)$/;
        $evalue = $1 if /^evalue=(.+)$/;
        $retries = $1 if /^retries=(\d+)$/;
        $maxHits = $1 if /^maxHits=(\d+)$/;
        $dbBlast = $1 if /^dbBlast=(.+)$/;
        $eqText = $1 if /^eqText=(.+)$/;
        $blastServer = $1 if /^blastServer=(.+)$/;
        $blastdbcmd = $1 if /^blastdbcmd=(.+)$/;
        $blastp = $1 if /^blastp=(.+)$/;
        $threads = $1 if /^threads=(\d+)$/;
    }
    close IN;
}


## read query list ##

open IN, "<$wkDir/result/bbh_input.txt" or die "Error: input file not accessible.\n";
while(<IN>) {
    s/\s+$//;
    @a = split (/\t/);
    # '1' is subject and '2' is query
    %h = ('set', $a[0], 'id', $a[1], 'accn1', $a[2], 'taxid1', $a[3], 'accn2', $a[4], 'taxid2', $a[5]);
    push @queries, {%h};
}
close IN;

## read previous results ##

if (-e "$wkDir/result/bbh.txt") {
    open IN, "<$wkDir/result/bbh.txt";
    while (<IN>) {
        s/\s+$//;
        @a = split (/\t/);
        $done{$a[0]."_".$a[2]} = 1;
    }
    close IN;
    my $nDone = 0;
    for ($i=0; $i<=$#queries; $i++) {
        if (exists $done{$queries[$i]{'set'}."_".$queries[$i]{'accn1'}}) {
            $queries[$i]{'done'} = 1;
            $nDone ++;
        }
    }
    if ((scalar @queries) <= $nDone) {
        print "Reversal BLAST is already completed. Exiting...\n";
        exit 0;
    } else {
        print "".(scalar @queries)." proteins in total, $nDone completed, remaining ".(scalar @queries - $nDone)." to BLAST.\n";
    }
} else {
    print "".(scalar @queries)." proteins to BLAST.\n";
}

## perform batch BLAST ##

for ($i=0; $i<=$#queries; $i++) {
    next if exists $queries[$i]{'done'};
    $isBBH = 0;
    $query = $queries[$i]{'accn2'};
    $query =~ s/\.[\d]+$//;
    $subject = $queries[$i]{'accn1'};
    $subject =~ s/\.[\d]+$//;
    $taxid = $queries[$i]{'taxid1'};
    $eqText = "txid$taxid [ORGN]";
    print "BLASTing $query ...";
    blast;
    open OUT, ">>$wkDir/result/bbh.txt";
    print OUT $queries[$i]{'set'}."\t".$queries[$i]{'id'}."\t".$queries[$i]{'accn1'}."\t".$queries[$i]{'taxid1'}."\t".$queries[$i]{'accn2'}."\t".$queries[$i]{'taxid2'}."\t$isBBH\n";    close OUT;
    print " $isBBH.\n";
}
unlink "$wkDir/blast.seq";
print "Reversal BLAST is completed.\n";

exit 0;


## subroutines ##

## perform BLAST ##

sub blast () {

    ## BLAST using standalone ncbi-blast+ program ##

    if ($httpBlast) {
        my @hits;
        # the report contains: accession, gi, length, taxid, sequence, title
        my @out = `$blastdbcmd -db $dbBlast -entry $query -outfmt \"%a %g %l %T %s %t\"`;
        foreach (@out) {
            s/\s+$//; @b = split (/\s+/);
            last if ($b[0] eq $query or $b[0] =~/^$query\.\d+$/);
        }
        my $length = $b[2];
        open TMP, ">$wkDir/blast.seq"; print TMP $b[4]; close TMP;
        # the report contains: subject accessions (all), E-value, bit score, aligned part of subject sequence
        @out = `$blastp -query $wkDir/blast.seq -db $dbBlast -entrez_query \"$eqText\" -remote -outfmt \"6 sallacc evalue bitscore\"`;
        @a = split (/\t/, $out[0]);
        my @accns = split (/;/, $a[0]);
        foreach (@accns) {
            if ($_ eq $subject) { $isBBH = 1; last; }
        }
    }

    ## BLAST via http connection to NCBI server ##

    else{
        $isError = 0;
        my $args = "CMD=Put&PROGRAM=blastp&DATABASE=$dbBlast&QUERY=$query&EQ_TEXT=$eqText";
        my $req = new HTTP::Request POST => $blastServer;
        $req->content_type('application/x-www-form-urlencoded');
        $req->content($args);
        my $response = $ua->request($req);
        my $rid;
        if ($response->content =~ /^    RID = (.*$)/m) { $rid = $1; } else { retry; return; };
        if ($response->content =~ /^    RTOE = (.*$)/m) { sleep $1; } else { retry; return; };
        while (1) {
            sleep 5;
            $args = "$blastServer?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
            $req = new HTTP::Request GET => $args;
            $response = $ua->request($req);
            if ($response->content =~ /\s+Status=WAITING/m) { next; }
            if ($response->content =~ /\s+Status=FAILED/m) { $isError = 1; last; }
            if ($response->content =~ /\s+Status=UNKNOWN/m) { $isError = 1; last; }
            if ($response->content =~ /\s+Status=READY/m) {
                if ($response->content =~ /\s+ThereAreHits=yes/m) { last; }
                else{ last; } # no hits;
            }
            $isError = 1;
            last;
        }
        if ($isError) { retry; return; }
        $args = "$blastServer?CMD=Get&ALIGNMENT_VIEW=Tabular&FORMAT_TYPE=Text&RID=$rid";
        $req = new HTTP::Request GET => $args;
        $response = $ua->request($req);
        if ($response->content !~ /blastp/) { retry; return; }
        my @out = split(/\n/, $response->content);
        my $read = 0;
        foreach (@out) {
            if (/hits? found/) { $read = 1; next; }
            next unless $read;
            @a = split (/\t/);
            if ($a[1] =~ $subject) { $isBBH = 1; }
            last;
        }
        unless ($read) { retry; return; }
    }
    return;
}

## retry BLAST ##

sub retry () {
    if ($iRetry < $retries) { # retry
        print ".";
        $iRetry ++;
        sleep 10;
        blast;
    } else { # fail
        $iRetry = 0;
        $isBBH = "error";
    }
}
