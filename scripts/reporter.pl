#!/usr/bin/perl

use warnings;
use strict;
$| = 1;

print "
-> Reporter: Generate reports of HGT prediction results. <-
";

print "

Usage:
  perl reporter.pl <working directory>

Output:
  result/detail/<name>.txt
  result/statistics/fingerprint.txt (optional)
  result/statistics/<name>.txt (optional)

" and exit unless @ARGV;


## all-purpose variables ##

my $i; my $j; my $n; my $s; my $t; my @a; my @b; my @c; my %h;


## program parameters ##

my $wkDir = $ARGV[0];
my $interactive = 1;

my $deHypo = 0;                    # ignore hypothetical proteins

my $byDonor = 1;                    # summarize HGT events by organism
my @ranks = ('species', 'genus', 'family', 'order', 'class', 'phylum');
my $sumRank = "order";            # on this rank
my $defOnly = 1;                    # ignore those without this rank defined

my $byFunction = 0;                # summarize HGT by function
my $dirFunction;                    # directory containing functional annotations

my $byOrthology = 0;                # generate report by ortholog
my $smOrthology;                    # file containing scheme of gene orthology
my $nameOGs = 1;                    # Name COGs if necessary
my $exORFan = 1;                    # Exclude ORFans
my @cogs;                            # clusters of orthologs (COGs) [ID] -> name short (long), predicted?
my %cogids = ();                    # clusters of orthologs (COGs): accn -> ID
my $isCOG = 0;                        # whether use COG
my %ORFans = ();                    # ORFans set -> (array)

my $outText = 1;                    # generate report in plain text
my $outHTML = 1;                    # generate report in web page (HTML)
my $outExcel = 0;                    # generate report in Excel spreadsheet
my $detailExcel = 1;                # attach detailed output in Excel workbook


## global variables ##

my $title;                            # title of the analysis
my @sets;                            # protein sets (genomes)

my %results = ();                    # master variable of results: $results{$set}{$accn}{$match}

my %totals = ();                    # total number of genes per genome
my %organisms = ();

my %taxadb = ();                    # taxa.db
my %ranksdb = ();                    # ranks.db

my $workbook;                        # Excel workbook
my $worksheet;                        # Excel worksheet
my $excelRow;                        # active row number
my $excelTitle;                    # Excel title format
my $excelHeader;                    # Excel header format
my ($excelGrey, $excelGreen, $excelYellow, $excelRed);                # Excel data formats

my $iLoss = 0; my $iPOE = 0; my $iMatch = 0;


## Read configuration ##

if (-e "$wkDir/config.txt") {
    open IN, "<$wkDir/config.txt";
    while (<IN>) {
        s/#.*$//; s/\s+$//g; s/^\s+//g; next unless $_;
        @sets = split (/,/, $1) if /^inSets=(.+)$/;
        $title = $1 if /^title=(.+)$/;
        $deHypo = $1 if /^deHypo=([012])$/;
        $detailExcel = $1 if /^detailExcel=([01])$/;

        $byDonor = $1 if /^byDonor=([01])$/;
        @ranks = split (/,/, $1) if /^ranks=(.+)$/;
        $sumRank = $1 if /^sumRank=(.+)$/;
        $defOnly = $1 if /^definedOnly=([01])$/;

        $byOrthology = $1 if /^byOrthology=([01])$/;
        $smOrthology = $1 if /^smOrthology=(.+)$/;
        $nameOGs = $1 if /^nameOGs=([01])$/;
        $exORFan = $1 if /^exORFan=([01])$/;

        $byFunction = $1 if /^byFunction=([01])$/;
        $dirFunction = $1 if /^dirFunction=(.+)$/;

        $outText = $1 if /^outText=([01])$/;
        $outExcel = $1 if /^outExcel=([01])$/;
        $outHTML = $1 if /^outHTML=([01])$/;

        $interactive = $1 if /^interactive=([01])$/;
    }
}


## check previous result ##

if (-e "$wkDir/result/summary.txt" or -d "$wkDir/result/HGT") {
    print "Warning: Summarized report from a previous analysis is detected.\n";
    if ($interactive) { print "Press Enter to overwrite, or press Ctrl+C to exit:\n"; $s = <STDIN>; }
    else{ print "To be overwritten.\n"; }
}


## Read input protein sets ##

unless (@sets) {
    opendir (DIR, "$wkDir/result/detail") or die "Error: result directory not accessible.\n";
    @sets = grep(/\.txt$/,readdir(DIR)) or die "Error: no result found.\n";
    for ($i=0; $i<=$#sets; $i++) { $sets[$i] =~ s/\.txt$//; }
    close DIR;
}


## Read taxonomy information ##

open IN, "<$wkDir/taxonomy/taxa.db";
while (<IN>) {
    s/\s+$//; next if /^#/; next unless $_;
    @a = split /\t/;
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


## Computing and/or reading COG scheme ##

if ($byOrthology) {
    if ($smOrthology) {
        unless (-e $smOrthology) {
            print "Orthology scheme file $smOrthology not accessible.\n";
            $smOrthology = "$wkDir/taxonomy/orthology.db";
        }
    } else { $smOrthology = "$wkDir/taxonomy/orthology.db"; }
    unless (-e $smOrthology) { # Identify orthology scenario using BBH
        $s = $0; $s =~ s/analyzer\.pl$/orthologer.pl/;
        system "$^X $s $wkDir";
    }
    unless (-e $smOrthology) {
        print "Error: Identification of orthology failed.\n";
        exit 1;
    }
    if ($nameOGs) { # Name OGs if necessary
        $i = 0;
        open IN, "<$smOrthology";
        while (<IN>) {
            s/\s+$//; next if /^#/; next unless $_;
            if (/\t/) { $i = 1; }
            last;
        }
        close IN;
        unless ($i) {
            system "perl scripts/cogNamer.pl $wkDir $smOrthology $wkDir/taxonomy/orthology.named.db";
        }
        $smOrthology = "$wkDir/taxonomy/orthology.named.db" if -e "$wkDir/taxonomy/orthology.named.db";
    }
    open IN, "<$smOrthology";
    while (<IN>) {
        s/\s+$//; next if /^#/; next unless $_;
        @a = split (/\t/);
        @b = split (/[\s\/,]/, $a[$#a]);
        next if ($exORFan and !$#b);
        %h = ();
        $h{'name'} = "";
        $h{'name'} = $a[0] if ($#a);
        $h{'data'} = [("") x @sets];
        push @cogs, {%h};
        $cogids{$_} = $#cogs for @b;
    }
    close IN;
    $isCOG = 1 if @cogs;
}


## Read prediction results ##

foreach my $set (@sets) {
    my %result = ();
    my %titles = ();
    @a = (); $ORFans{$set} = [@a];
    open IN, "<$wkDir/result/detail/$set.txt";
    while (<IN>) {
        s/\s+$//; next unless $_; next if (/HGTector/);
        my @a = split (/\t/);
        if (/^Query\t/) {
            next if ($iLoss+$iPOE+$iMatch);
            for ($i=0; $i<=$#a; $i++) {
                $iLoss = $i if ($a[$i] eq "Loss");
                $iPOE = $i if ($a[$i] eq "POE");
                $iMatch = $i if ($a[$i] eq "Match");
            }
            next;
        }
        if ($deHypo) { # ignore hypothetical proteins
            if ($isCOG and exists $cogids{$a[0]}) {
                my $name = $cogs[$cogids{$a[0]}]{'name'};
                next if ($name =~ /hypothetical/ or $name =~ /hypotethical/ or $name =~ /hypothetcial/);
            } else {
                next if ($a[2] =~ /hypothetical/ or $a[2] =~ /hypotethical/ or $a[2] =~ /hypothetcial/);
            }
        }
        my %gene = ();
        $gene{'length'} = $a[1];
        $gene{'product'} = $a[2];
        $gene{'hits'} = $a[3];
        if ($isCOG) {
            if (exists $cogids{$a[0]}) {
                $gene{'cogid'} = $cogids{$a[0]};
            } else {
                push @{$ORFans{$set}}, "$a[0]|$a[2]";
            }
        }
        if ($#a > 4) {
            $gene{'self'} = $a[4];
            $gene{'close'} = $a[5];
            $gene{'distal'} = $a[6];
            $gene{'hgt'} = $a[7];
            $gene{'loss'} = $a[$iLoss] if $iLoss;
            $gene{'poe'} = $a[$iPOE] if $iPOE;
            $gene{'match'} = $a[$iMatch] if $iMatch;
        }
        $result{$a[0]} = {%gene};
    }
    close IN;
    $results{$set} = {%result};
    @a = keys %result;
    $totals{$set} = @a;
}


## start to generate reports ##

@a = localtime(time);
my $now = (++$a[4])."/$a[3]/".($a[5]+=1900)." $a[2]:$a[1]:$a[0]";
$title = "HGTector report at $now" unless $title;

## load module to access Excel ##

if ($outExcel) {
    eval { require Spreadsheet::WriteExcel; Spreadsheet::WriteExcel->import() };
    print "Error: Perl module Spreadsheet::WriteExcel is not installed.\n" and exit 1 if ($@);
    $workbook  = Spreadsheet::WriteExcel->new("$wkDir/result/report.xls");
    $excelTitle = $workbook->add_format (bold=>1, size=>12, valign=>'vcenter');
    $excelHeader = $workbook->add_format(bold=>1, align=>'center', top=>1, bottom=>1);
    $excelGrey = $workbook->add_format(align=>'center', fg_color=>'silver', pattern=>2, border=>1, border_color=>'white');
    $excelGreen = $workbook->add_format(align=>'center', fg_color=>'lime', pattern=>4, border=>1, border_color=>'white');
    $excelYellow = $workbook->add_format(align=>'center', fg_color=>'yellow', pattern=>2, border=>1, border_color=>'white');
    $excelRed = $workbook->add_format(align=>'center', fg_color=>'red', pattern=>4, border=>1, border_color=>'white');
## This is a good example of conditional formatting code, but I cannot use more than three rules.
#    $excelData = $workbook->add_format (num_format=>'[Red][<=0]0;[Green][<=20]0;[Blue]0');
}


if ($outHTML) {
    open HTML, ">$wkDir/result/report.html";
    print HTML "<!DOCTYPE html>\n<html>\n <head>\n  <title>$title</title>\n </head>\n <body>\n";
    print HTML "  <p><h1>$title</h1></p>\n";
}


## Generate general report ##

mkdir "$wkDir/result/HGT" unless -d "$wkDir/result/HGT";
if ($iLoss) { mkdir "$wkDir/result/loss" unless -d "$wkDir/result/loss"; }
if ($iPOE) { mkdir "$wkDir/result/POE" unless -d "$wkDir/result/POE"; }
if ($outText) {
    open OUT, ">$wkDir/result/summary.txt";
}
if ($outHTML) {
    print HTML "  <p><h2>Summary</h2></p>";
    print HTML "  <table border=1>\n   <tr>\n    <th width=120>Genome</th>\n    <th width=120>Total</th>\n    <th width=120>HGT</th>\n";
    print HTML "    <th width=120>Loss</th>\n" if $iLoss;
    print HTML "    <th width=120>POE</th>\n" if $iPOE;
    print HTML "   </tr>\n";
}
if ($outExcel) {
    $worksheet = $workbook->add_worksheet ('Summary');
    $worksheet->set_row (0, 24);
    $worksheet->set_column (0, 4, 12);
    $worksheet->write (0, 0,  $title, $excelTitle);
    $worksheet->write (1, 0, ['Genome', 'Total', 'HGT'], $excelHeader);
    $worksheet->write (1, 3, 'Loss', $excelHeader) if $iLoss;
    $i = 3; $i ++ if $iLoss;
    $worksheet->write (1, $i, 'POE', $excelHeader) if $iPOE;
    $excelRow = 2;
}

my $allhgt; # total number of HGT-derived genes
foreach my $set (@sets) {
    my $hgt = 0; my $loss = 0; my $poe = 0;
    open HGT, ">$wkDir/result/HGT/$set.txt";
    open LOSS, ">$wkDir/result/loss/$set.txt" if $iLoss;
    open POE, ">$wkDir/result/POE/$set.txt" if $iPOE;
    my %result = %{$results{$set}};
    my $total = keys %result;
    foreach my $accn (sort keys %result) {
        print HGT $accn."\n" and $hgt ++ if (exists $result{$accn}{'hgt'} and $result{$accn}{'hgt'});
        if ($iLoss) { print LOSS $accn."\n" and $loss ++ if (exists $result{$accn}{'loss'} and $result{$accn}{'loss'}); }
        if ($iPOE) { print POE $accn."\n" and $poe ++ if (exists $result{$accn}{'poe'} and $result{$accn}{'poe'}); }
    }
    close HGT; close LOSS if $iLoss; close POE if $iPOE;
    if ($outText) {
        print OUT "$set has $total "; print OUT "non-hypothetical " if $deHypo; print OUT "protein-coding genes.\n";
        print OUT "  HGT: $hgt";
        print OUT ", Loss: $loss" if $iLoss;
        print OUT ", POE: $poe" if $iPOE;
        print OUT ".\n\n";
    }
    if ($outHTML) {
        print HTML "   <tr>\n    <td>$set</td>\n    <td><a href=\"detail/$set.txt\">$total</a></td>\n    <td><a href=\"HGT/$set.txt\">$hgt</a></td>\n";
        print HTML "    <td><a href=\"loss/$set.txt\">$loss</a></td>\n" if $iLoss;
        print HTML "    <td><a href=\"POE/$set.txt\">$poe</a></td>\n" if $iPOE;
        print HTML "   </tr>\n";
    }
    if ($outExcel) {
        $worksheet->write ($excelRow, 0, [$set, $total, $hgt]);
        $worksheet->write ($excelRow, 3, $loss) if $iLoss;
        $i = 3; $i ++ if $iLoss;
        $worksheet->write ($excelRow, $i, $poe) if $iPOE;
        $excelRow ++;
    }
    $allhgt += $hgt;
}
close OUT if ($outText);
print HTML "  </table>\n" if ($outHTML);


## Generate by donor group report ##

if ($byDonor) {
    foreach my $set (@sets) {
        my %result = %{$results{$set}};
        my $count = 0;
        foreach my $accn (keys %result) {
            $count ++;
            next unless exists $result{$accn}{'hgt'};
            next unless $result{$accn}{'hgt'};
            next unless exists $result{$accn}{'match'};
            next unless $result{$accn}{'match'};
            my $organism = $result{$accn}{'match'};
            my $group = "";
            if ($sumRank) {
                $organism =~ /^(\d+) \(/;
                if (exists $taxadb{$1} and exists $taxadb{$1}{$sumRank} and $taxadb{$1}{$sumRank} and exists $ranksdb{$taxadb{$1}{$sumRank}}) {
                    # next unless exists $taxadb{$1}{'class'} and $taxadb{$1}{'class'} eq '28211'; # alphaproteobacteria
                    # next unless substr ($taxadb{$1}{'rank'}, 0, 3) eq "/2/"; # bacteria
                    # next if exists $taxadb{$1}{'class'} and $taxadb{$1}{'class'};
                    # $organism = $ranksdb{$taxadb{$1}{'phylum'}}.",".$ranksdb{$taxadb{$1}{'class'}}.",".$ranksdb{$taxadb{$1}{'order'}};
                    $organism = $ranksdb{$taxadb{$1}{$sumRank}};
                } else { next; }
                # if (exists $taxadb{$1} and substr ($taxadb{$1}{'rank'}, 0, 3) eq "/2/") { $organism = "Bacteria"; }
                # elsif (exists $taxadb{$1} and substr ($taxadb{$1}{'rank'}, 0, 6) eq "/2157/") { $organism = "Archaea"; }
                # elsif (exists $taxadb{$1} and substr ($taxadb{$1}{'rank'}, 0, 6) eq "/2759/") { $organism = "Eukaryota"; }
                # else{ next; }
            }
            $results{$set}{$accn}{'group'} = $organism;
            if (exists $organisms{$organism}) {
                if ($organisms{$organism}{$set}) { $organisms{$organism}{$set} .= ",$accn"; }
                else { $organisms{$organism}{$set} .= $accn; }
            } else {
                my %h = ();
                $h{$_} = "" for (@sets);
                $h{$set} = $accn;
                $organisms{$organism} = {%h};
            }
        }
    }

    # head row
    if ($outText) {
        open OUT, ">$wkDir/result/donor.txt";
        print OUT "Organism\tMean\t".join ("\t", @sets)."\n";
    }
    if ($outHTML) {
        print HTML "  <br><hr>\n  <p></p>\n  <p><h2>HGT by putative donor group</h2></p>\n";
        print HTML "  <p><font size='-1'>HGT-derived genes are summarized by putative donor group as indicated by the best distal match. Note that the real donor might be an ancestor of the match organism, therefore, a higher taxonomic rank is recommended to describe the donor group.";
        print HTML "<br><font style='background-color: lightgreen'>green</font>: 1-4, <font style='background-color: khaki'>yellow</font>: 5-9, <font style='background-color: lightpink'>red</font>: 10+</font></p>\n";
        print HTML "  <table border=1>\n   <tr>\n    <th width=250>Donor</th>\n    <th width=80>Mean</th>\n";
        print HTML "    <th width=120>$_</th>\n" for (@sets);
        print HTML "   </tr>\n";
    }
    if ($outExcel) {
        $worksheet = $workbook->add_worksheet ('Donor');
        $worksheet->set_row (0, 24);
        $worksheet->set_column (0, 0, 16);
        $worksheet->set_column (1, scalar @sets, 7);
        $worksheet->write (0, 0,  "HGT by donor organism", $excelTitle);
        $worksheet->write (1, 0, "Donor", $excelHeader);
        $worksheet->write (1, 1, "Mean", $excelHeader);
        $worksheet->write (1, 2, \@sets, $excelHeader);
        $excelRow = 3;
    }

    # total row
    @a = (); push (@a, $totals{$_}) for (@sets);
    $j = 0; $j += $_ for (@a); $j = sprintf("%d", $j/@a);
    if ($outText) {
        print OUT "Total\t$j\t".join ("\t", @a)."\n";
    }
    if ($outHTML) {
        print HTML "   <tr>\n    <td>Total</td>\n";
        print HTML "    <td>$j</td>\n";
        print HTML "    <td>$_</td>\n" for (@a);
        print HTML "   </tr>\n";
    }
    if ($outExcel) {
        $worksheet->write (2, 0, "Total");
        $worksheet->write (2, 1, $j);
        $worksheet->write (2, 2, \@a);
    }

    # table content
    foreach my $organism (sort keys %organisms) {
        my @accns = (); # protein accns (comma-separated)
        my @nums = (); # numbers only
        foreach (@sets) {
            $s = $organisms{$organism}{$_};
            push (@accns, $s);
            @a = split (/,/, $s);
            $s = @a; push (@nums, $s);
        }
        $j = 0; $j += $_ for (@nums); $j = sprintf("%d", $j/@nums);
        if ($outText) {
            print OUT "$organism\t$j\t".join("\t", @nums)."\n";
        }
        if ($outHTML) {
            print HTML "   <tr>\n    <td>$organism</td>\n    <td>$j</td>\n";
            for ($i=0; $i<=$#nums; $i++) {
                $n = $nums[$i];
                $t = "";
                $s = join ("&#013;", split (/,/, $accns[$i]));
                if ($n > 0 and $n < 5) { $t = "lightgreen"; }
                elsif ($n >= 5 and $n < 10) { $t = "khaki"; }
                elsif ($n >= 10) { $t = "lightpink"; }
                print HTML "    <td title='$s' style='background-color: $t;'>$n</td>\n";
            }
            print HTML "   </tr>\n";
        }
        if ($outExcel) {
            $worksheet->write ($excelRow, 0, $organism);
            $worksheet->write ($excelRow, 1, $j);
            for ($i=0; $i<=$#nums; $i++) {
                next unless $nums[$i]; # (=0)
                if ($nums[$i] < 5) { $worksheet->write_number ($excelRow, $i+2, $nums[$i], $excelGreen); }
                elsif ($nums[$i] < 10) { $worksheet->write_number ($excelRow, $i+2, $nums[$i], $excelYellow); }
                else{ $worksheet->write_number ($excelRow, $i+2, $nums[$i], $excelRed); }
                $worksheet->write_comment ($excelRow, $i+2, join ("\n", split (/,/, $accns[$i])));
            }
            $excelRow ++;
        }
    }
    close OUT if ($outText);
    print HTML "  </table>\n" if ($outHTML);
    print "Report by donor organism generated.\n";
}


## Generate by function report ##

if ($byFunction and -d $dirFunction) {

    ## read Blast2GO output ##

    my %goes = ();                            # master record of GOes. GOes{ID} = { $group, $term, $nGene, $nHGT, $genes}
    opendir DIR, $dirFunction;
    my @files = readdir(DIR);
    close DIR;
    foreach my $file (@files) {
        my $set = "";
        if ($file =~ /(.+)\.[^.]+$/) { $set = $1; }
        else{ $set = $file; }
        next unless exists $results{$set};
        open IN, "<$dirFunction/$file";
        while (<IN>) {
            s/\s+$//;
            next if /^SeqName/;
            @a = split (/\t/);
            push @a, "" if $#a == 3;
            my $gene = "";
            if ($a[0] =~ /\|/) {
                @b = split (/\|/, $a[0]);
                $b[$#b] =~ s/\.\d+$//;
                if ($b[$#b] =~ /cdsid_(.+)$/) { $gene = $1; }
                else{ $gene = $b[$#b]; }
            } else {
                $a[0] =~ s/\.\d+$//;
                $gene = $a[0];
            }
            next unless $gene;
            next unless exists $results{$set}{$gene};
            $a[3] =~ s/^GO://;
            if (exists $goes{$a[3]}) {
                $goes{$a[3]}{$set}{'genes'} .= ",$gene";
            } else {
                %h = ('nGene', 0, 'nHGT', 0, 'genes', "");
                my %go = ('group', $a[2], 'term', $a[4]);
                $go{$_} = {%h} for (@sets);
                $go{$set}{'genes'} = $gene;
                $goes{$a[3]} = {%go};
            }
            $goes{$a[3]}{$set}{'nGene'} ++;
            $goes{$a[3]}{$set}{'nHGT'} ++ if $results{$set}{$gene}{'hgt'};
            if (exists $results{$set}{$gene}{'goes'}) { $results{$set}{$gene}{'goes'} .= "; $a[3] ($a[2]) $a[4]"; }
            else{ $results{$set}{$gene}{'goes'} = "$a[3] ($a[2]) $a[4]"; }
        }
        close IN;
    }

    ## write report by function ##

    if ($outText) {
        open OUT, ">$wkDir/result/function.txt";
        print OUT "GO\tGroup\tTerm".join ("\t", @sets)."\n";
    }
    if ($outHTML) {
        print HTML "  <br><hr>\n  <p></p>\n  <p><h2>HGT by functional annotation</h2></p>";
        print HTML "  <p><font size='-1'>HGT-derived genes are summarized by annotated gene ontology (GO). Cell values are ratio of HGT-derived genes associated with a GO versus all genes associated with this GO. Every GO is counted if multiple GOes are associated with one gene.";
        print HTML "<br>white: 0, <font style='background-color: lightgreen'>green</font>: (0, 0.1), <font style='background-color: khaki'>yellow</font>: [0.1, 0.5), <font style='background-color: lightpink'>red</font>: [0.5, 1]</font></p>\n";
        print HTML "  <table border=1>\n   <tr>\n    <th width=60>GO</th>\n    <th width=240>Term</th>\n";
        print HTML "    <th width=120>$_</th>\n" for (@sets);
        print HTML "   </tr>\n";
    }
    if ($outExcel) {
        $worksheet = $workbook->add_worksheet ('Function');
        $worksheet->set_row (0, 24);
        $worksheet->set_column (0, 0, 8);
        $worksheet->set_column (1, 1, 40);
        $worksheet->set_column (2, 1 + scalar @sets, 7);
        $worksheet->write (0, 0,  "HGT by function", $excelTitle);
        $worksheet->write (1, 0, ["GO","Term"], $excelHeader);
        $worksheet->write (1, 2, \@sets, $excelHeader);
        $excelRow = 2;
    }
    foreach my $group ('C', 'P', 'F') {
        if ($outHTML) {
            if ($group eq "C") { print HTML "    <tr><th colspan = '2'>Cellular component</th></tr>\n"; }
            if ($group eq "P") { print HTML "    <tr><th colspan = '2'>Biological process</th></tr>\n"; }
            if ($group eq "F") { print HTML "    <tr><th colspan = '2'>Molecular function</th></tr>\n"; }
        }
        if ($outExcel) {
            if ($group eq "C") { $worksheet->write ($excelRow++, 1, "Cellular component", $excelHeader); }
            if ($group eq "P") { $worksheet->write ($excelRow++, 1, "Biological process", $excelHeader); }
            if ($group eq "F") { $worksheet->write ($excelRow++, 1, "Molecular function", $excelHeader); }
        }
        foreach my $go (sort {$goes{$a}{'term'} cmp $goes{$b}{'term'}} keys %goes) {
            next unless $goes{$go}{'group'} eq $group;
            next if (($goes{$go}{'term'} eq "cellular_component") or ($goes{$go}{'term'} eq "biological_process") or ($goes{$go}{'term'} eq "molecular_function"));
            if ($outText) { print OUT "$go\t$goes{$go}{'group'}\t$goes{$go}{'term'}"; }
            if ($outHTML) { print HTML "   <tr>\n    <td>$go</td>\n    <td>$goes{$go}{'term'}</td>\n"; }
            if ($outExcel) {
                $worksheet->write_string ($excelRow, 0, $go);
                $worksheet->write_string ($excelRow, 1, $goes{$go}{'term'});
            }
            $j = 2;
            foreach (@sets) {
                if ($goes{$go}{$_}{'nGene'} and $goes{$go}{$_}{'nHGT'}) { $i = sprintf ("%.2g", $goes{$go}{$_}{'nHGT'}/$goes{$go}{$_}{'nGene'}); }
                else { $i = 0; }
                if ($outText) { print OUT "\t$goes{$go}{$_}{'nHGT'}\/$goes{$go}{$_}{'nGene'}"; }
                if ($outHTML) {
                    print HTML "    <td title='$goes{$go}{$_}{'nHGT'}\/$goes{$go}{$_}{'nGene'}'";
                    if ($i) {
                        if ($i > 0 and $i < 0.1) { $t = "lightgreen"; }
                        elsif ($i >= 0.1 and $i < 0.5) { $t = "khaki"; }
                        elsif ($i >= 0.5) { $t = "lightpink"; }
                        print HTML " style='background-color: $t;'>$i</td>\n";
                    } else {
                        print HTML ">$i</td>\n";
                    }
                }
                if ($outExcel) { $worksheet->write_number ($excelRow, $j, $i); $worksheet->write_comment ($excelRow, $j++, "$goes{$go}{$_}{'nHGT'}\/$goes{$go}{$_}{'nGene'}"); }
            }
            if ($outText) { print OUT "\n"; }
            if ($outHTML) { print HTML "   </tr>\n"; }
            if ($outExcel) { $excelRow++; }
        }
    }
    if ($outText) { close OUT; }
    if ($outHTML) { print HTML "  </table>\n";}
    print "Report by functional annotation generated.\n";
}


## Generate by ortholog report ##

if ($byOrthology) {
    for ($i=0; $i <=$#sets; $i++) {
        foreach my $accn (keys %{$results{$sets[$i]}}) {
            %h = %{$results{$sets[$i]}{$accn}};
            if (exists $h{'cogid'}) {
                if (exists $h{'hgt'} and $h{'hgt'}) {
                    $cogs[$h{'cogid'}]{'data'}[$i] .= "1:$accn ";
                } else {
                    $cogs[$h{'cogid'}]{'data'}[$i] .= "0:$accn ";
                }
            }
        }
    }
    if ($outText) {
        open OUT, ">$wkDir/result/orthology.txt";
        print OUT "ID\tName\t".join ("\t", @sets)."\n";
    }
    if ($outHTML) {
        print HTML "  <br><hr>\n  <p></p>\n  <p><h2>HGT by gene orthology</h2></p>";
        print HTML "  <p><font size='-1'>HGT-derived genes are summarized by orthologous groups (OGs).";
        print HTML "<br>empty: gene not present, <font style='background-color: lightgrey'>0</font>: gene present but not predicted as HGT, <font style='background-color: lightgreen'>1</font>: gene present and predicted as HGT, <font style='background-color: khaki'>m/n</font> (n>=2): n paralogs are present, in which m are predicted as HGT.</font></p>\n";
        if ($cogs[0]{'name'}) { $i = 250; } else { $i = 50; }
        print HTML "  <table border=1>\n   <tr>\n    <th width=$i>OG</th>\n";
        print HTML "    <th width=120>$_</th>\n" for (@sets);
        print HTML "   </tr>\n";
    }
    if ($outExcel) {
        $worksheet = $workbook->add_worksheet ('Orthology');
        $worksheet->set_row (0, 24);
        if ($cogs[0]{'name'}) { $i = 40; } else { $i = 6; }
        $worksheet->set_column (0, 0, $i);
        $worksheet->set_column (1, scalar @sets, 7);
        $worksheet->write (0, 0,  "HGT by gene orthology", $excelTitle);
        $worksheet->write (1, 0, "OG", $excelHeader);
        $worksheet->write (1, 1, \@sets, $excelHeader);
        $excelRow = 2;
    }

    my $allCOG; my $hgtCOG;

    for ($i=0; $i <=$#cogs; $i++) {
        @a = @{$cogs[$i]{'data'}};
        $s = 0; foreach (@a) { if (/^1/) { $s = 1; last; } } next unless $s;

        for ($j=0; $j<=$#a; $j++) {
            next unless $a[$j];
            @b = split (/\s/, $a[$j]);
            $allCOG += scalar @b;
            foreach (@b) { $hgtCOG ++ if /^1/; }
        }

        if ($outText) {
            print OUT "$i\t$cogs[$i]{'name'}\t".join("\t",@a)."\n";
        }
        if ($outHTML) {
            if ($cogs[$i]{'name'}) { $s = $cogs[$i]{'name'}; } else { $s = $i; }
            print HTML "   <tr>\n    <td>$s</td>\n";
            for ($j=0; $j<=$#a; $j++) {
                if ($a[$j]) {
                    @b = split (/\s/, $a[$j]);
                    my $hgt = 0; my $first; $t = "";
                    foreach (@b) {
                        /^([01]):(.+)$/;
                        $first = $2 unless $first;
                        $t .= "$2,";
                        $hgt ++ if $1 == 1;
                    }
                    print HTML "    <td style='background-color: ";
                    if ($#b) { print HTML "khaki"; }
                    elsif ($hgt) { print HTML "lightgreen"; }
                    else { print HTML "lightgrey"; }
                    print HTML ";' title='".join ("&#013;", split (/,/, $t))."'>";
                    print HTML "<a href='../search/$sets[$j]/$first.txt' >";
                    unless ($#b) { print HTML "$hgt</td>\n"; }
                    else{ print HTML "$hgt/".scalar(@b)."</td>\n"; }
                } else {
                    print HTML "    <td></td>\n";
                }
            }
            print HTML "   </tr>\n";
        }
        if ($outExcel) {
            if ($cogs[$i]{'name'}) { $s = $cogs[$i]{'name'}; } else { $s = $i; }
            $worksheet->write ($excelRow, 0, $s);
            for ($j=0; $j<=$#a; $j++) {
                next unless $a[$j];
                @b = split (/\s/, $a[$j]);
                my $hgt = 0; my $first; $t = "";
                foreach (@b) {
                    /^([01]):(.+)$/;
                    $first = $2 unless $first;
                    $t .= "$2,";
                    $hgt ++ if $1 == 1;
                }
                if ($#b) {
                    if ($hgt) { $worksheet->write_url ($excelRow, $j+1, $2, "$hgt/".scalar(@b), $excelGreen); }
                    else { $worksheet->write_url ($excelRow, $j+1, $2, "$hgt/".scalar(@b), $excelGrey); }
                } else {
                    if ($hgt) { $worksheet->write_url ($excelRow, $j+1, $2, 1, $excelGreen); }
                    else { $worksheet->write_url ($excelRow, $j+1, $2, 0, $excelGrey); }
                }
                $worksheet->write_comment ($excelRow, $j+1, $t);
            }
            $excelRow ++;
        }
    }
    print HTML "  </table>\n" if ($outHTML);

    if (%ORFans and 0) {
        if ($outText) { print OUT "#Singleton ORFans:\n"; }
        if ($outHTML) { print HTML "  <br><hr>\n  <p></p>\n  <p><h3>Singleton ORFans</h2></p>\n"; }
        foreach my $set (@sets) {
            if ($outText) { print OUT "$set: "; }
            if ($outHTML) { print HTML "    <p><b>$set</b>: "; }
            foreach my $gene (@{$ORFans{$set}}) {
                if ($outText) { print OUT "$gene "; }
                if ($outHTML) { print HTML "$gene "; }
            }
            if ($outText) { print OUT "\n"; }
            if ($outHTML) { print HTML "</p>\n"; }
        }
    }

    close OUT if ($outText);

    print "Report by gene orthology generated.\n";

    print "  Positive rate = $allhgt / ".($allhgt-$hgtCOG+$allCOG)." = ". (sprintf("%.3f", $allhgt/($allhgt-$hgtCOG+$allCOG))).".\n";

}



## Attach individual reports in spreadsheet ##

if ($outExcel and $detailExcel) {
    foreach my $set (@sets) {
        $worksheet = $workbook->add_worksheet ($set);
        $worksheet->set_row (0, 24);
        $worksheet->set_column (0, 0, 14);
        $worksheet->set_column (1, 1, 8);
        $worksheet->set_column (2, 2, 32);
        $worksheet->set_column (3, 9, 8);
        $worksheet->set_column ($iMatch, $iMatch+$byDonor, 32) if $iMatch;
        $worksheet->write (0, 0, $set, $excelTitle);
        $worksheet->write (1, 0, ["Query","Length","Product","Hits","Self","Close","Distal","HGT"], $excelHeader);
        $worksheet->write (1, $iLoss, "Loss", $excelHeader) if $iLoss;
        $worksheet->write (1, $iPOE, "POE", $excelHeader) if $iPOE;
        $worksheet->write (1, $iMatch, "Best distal match", $excelHeader) if $iMatch;
        $worksheet->write (1, $iMatch+1, "Putative donor group", $excelHeader) if $byDonor;
        $worksheet->write (1, $iMatch+$byDonor+1, "Function", $excelHeader) if $byFunction;
        $excelRow = 1;
        %h = %{$results{$set}};
        foreach my $accn (sort keys %h) {
            $excelRow ++;
            $worksheet->write_url ($excelRow, 0, "../search/$set/$accn.txt", $accn);
            $worksheet->write_number ($excelRow, 1, $h{$accn}{'length'});
            $worksheet->write ($excelRow, 2, $h{$accn}{'product'});
            $worksheet->write_number ($excelRow, 3, $h{$accn}{'hits'});
            next unless exists $h{$accn}{'self'};
            $worksheet->write_number ($excelRow, 4, $h{$accn}{'self'});
            $worksheet->write_number ($excelRow, 5, $h{$accn}{'close'});
            $worksheet->write_number ($excelRow, 6, $h{$accn}{'distal'});
            $worksheet->write_number ($excelRow, 7, $h{$accn}{'hgt'}) if $h{$accn}{'hgt'};
            $worksheet->write_number ($excelRow, $iLoss, $h{$accn}{'loss'}) if ($iLoss and $h{$accn}{'loss'});
            $worksheet->write_number ($excelRow, $iPOE, $h{$accn}{'poe'}) if ($iPOE and $h{$accn}{'poe'});
            next unless ($iMatch and $h{$accn}{'match'});
            $h{$accn}{'match'} =~ /^(\d+) \((.+)\)$/;
            $worksheet->write ($excelRow, $iMatch, $2);
            if ($byDonor and exists $h{$accn}{'group'}) {
                $worksheet->write ($excelRow, $iMatch+1, $h{$accn}{'group'});
            }
            if ($byFunction and exists $h{$accn}{'goes'}) {
                $s = $h{$accn}{'goes'};
                $s =~ s/GO:\d+ \(.\) //g;
                $worksheet->write ($excelRow, $iMatch+$byDonor+1, $s);
            }
            # $worksheet->write_comment ($excelRow, 10, $1);
        }
    }
}

if ($outHTML) {
    print HTML "\t</body>\n</html>\n";
    close HTML;
}
if ($outExcel) {
    $worksheet = $workbook->sheets(0);
    $worksheet->set_first_sheet();
    $worksheet->activate();
    $worksheet->set_selection(0, 0);
    $workbook->close();
}

exit 0;
