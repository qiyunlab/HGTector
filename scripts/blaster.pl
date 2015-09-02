#!/usr/bin/perl

use strict;
use warnings;
use threads;
use threads::shared;
$| = 1;


print "

This script performs BLAST in a batch mode with all the proteins given
in the protein list against the NCBI nr database.

Usage:
  Place input files in input/, one file per protein set, then execute:
  perl blaster.pl <working directory> <input file>

Types of input file:
  1. Proteins that are already recorded in GenBank.
    1a. NCBI protein summary file containing multiple protein records.
    1b. NCBI annotated genome file in GenBank format (*.gb or *.gbk).
    1c. A list of NCBI accession numbers (one entry per line).
  2. User-defined protein sequences in multifasta format.
    The sequence names should be unique and not containing the NCBI
    delimiter \"|\".

Output:
  blast/<name>/<accn>.bla
  blast/<name>.log

" and exit unless @ARGV;

print "
-> Blaster: Perform batch BLAST searches. <-
";


## public variables ##

my $i; my $j; my $n; my $s; my $t; my @a; my @b; my @c; my %h;

## modules ##

use LWP::Simple;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);
# use URI::Escape;

## global variables ##

my $wkDir = $ARGV[0];							# working directory
my $interactive = 1;							# interactive or automatic mode

my $nPt = 0;									# total number of proteins
my $nDone = 0;									# number of proteins that have beeb blasted

my $isInSeq = 0;								# use query sequences instead of accn
my %inSeqs = ();								# query protein sequences (if appicable)

my %ins = ();									# all input proteins

my $ua = LWP::UserAgent->new;

my %taxdumps = ();								# local taxonomy database: taxid => (name, parent, rank)
my %gi2taxids = ();								# GI to TaxID dictionary
my %dbTaxa :shared = ();						# taxa.db
my %dbRanks :shared = ();						# ranks.db
my %badTaxids = ();								# TaxIDs that don't exist in the local taxonomy database

## subroutines ##

sub http_blast;									# parameters: query (accn. no.), set; return: 1 - succeeded, 0 - fail
sub local_blast;								# same as above
sub self_blast;									# blast against itself
sub get_taxonomy;								# paramter: array of TaxIDs
sub order_accns;								# reorder accession number

## program parameters ##

my $blastMode = 0;								# BLAST mode (0: via http connection, 1: local BLAST)
my $preBlast = "";								# Folder containing pre-computed BLAST results.
my $selfTax = "";

my $nRequests = 1;								# multiple requests for http BLAST (default: 1)
my $nThreads = 1;								# multiple threads for local BLAST program
my $nQueries = 1;								# multiple queries per run for local BLAST program (0 for all)

my $nRetry = 5;									# maximum number of retries (http BLAST)
my $delay = 5;									# time delay (seconds) between two http requests (http BLAST)

my $nHits = 500;								# number of hits to return
my $maxHits = 0;								# maximum number of valid hits to preserve. if 0 then = nHits
my $evalue = 0.00001;							# maximum E-value cutoff
my $identity = 0;								# minimum percent identity cutoff
my $coverage = 0;								# minimum query coverage cutoff

my $dbBlast = "nr";
my $exUncultured = 1;							# exclude uncultured and environmental samples
my @searchTaxids = ();							# search under the following taxon groups (taxids)
my @ignoreTaxids = ();							# ignore organisms under the following taxids
my $eqText = "";								# entrez query parameter

my $seqBlast = 0;								# retrieve hit sequences
my $taxBlast = 0;								# retrieve taxonomy report
my $alnBlast = 0;								# retrieve multiple sequence alignment (conflicts seqBlast)

my $mergeDuplicates = 1;						# ignore hits with same taxon names and bit scores
my $taxonUCASE = 1;								# ignore taxon names that do not start with a capital letter
my @ignoreTaxa = ();							# ignore taxon names containing the following words
my $ignoreParalogs = 1;							# ignore potential paralogs (hits with same taxon names but different bit scores)
my $ignoreSeqRepeats = 1;						# ignore repeated sequences (hits targetting different regions of the same protein)
my $ignoreSubspecies = 0;						# ignore more than one subspecies from the same species
my @ranks = ('species', 'genus', 'family', 'order', 'class', 'phylum');

my $blastdbcmd = "blastdbcmd";
my $blastp = "blastp";
my $taxdump = "";								# directory of NCBI taxonomy database (nodes.dmp and names.dmp)
my $gi2taxid = "";								# GI to TaxID dictionary file

my $blastServer = "http://blast.ncbi.nlm.nih.gov/Blast.cgi";
my $blastAlignServer = "http://blast.ncbi.nlm.nih.gov/BlastAlign.cgi";
my $eSearchServer = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
my $eFetchServer = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
my $eSummaryServer = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi";

my $retryFailed = 1;


## read configurations ##

if (-e "$wkDir/config.txt"){
	open IN, "<$wkDir/config.txt";
	while (<IN>){
		s/#.*$//; s/\s+$//g; s/^\s+//g; next unless $_;
		$interactive = $1 if /^interactive=([01])$/;
		$blastMode = $1 if /^blastMode=(\d)$/;
		$preBlast = $1 if /^preBlast=(.+)$/;
		$preBlast =~ s/\/$//;
		$selfTax = $1 if /^selfTax=(.+)$/;
		
		$nHits = $1 if /^nHits=(\d+)$/;
		$evalue = $1 if /^evalue=(.+)$/;
		$identity = $1 if /^identity=(.+)$/;
		$coverage = $1 if /^coverage=(.+)$/;
		$identity = $1 if /^percIdent=(.+)$/; # backward compatibility

		$nRetry = $1 if /^nRetry=(\d+)$/;
		$delay = $1 if /^delay=(\d+)$/;
		$maxHits = $1 if /^maxHits=(\d+)$/;
		$dbBlast = $1 if /^dbBlast=(.+)$/;
		$eqText = $1 if /^eqText=(.+)$/;
		
		$blastServer = $1 if /^blastServer=(.+)$/;
		$blastAlignServer = $1 if /^blastAlignServer=(.+)$/;
		$eSearchServer = $1 if /^eSearchServer=(.+)$/;
		$eFetchServer = $1 if /^eFetchServer=(.+)$/;
		$eSummaryServer = $1 if /^eSummaryServer=(.+)$/;
		
		$blastdbcmd = $1 if /^blastdbcmd=(.+)$/;
		$blastp = $1 if /^blastp=(.+)$/;
		$taxdump = $1 if /^taxdump=(.+)$/;
		$taxdump =~ s/\/$//;
		$gi2taxid = $1 if /^gi2taxid=(.+)$/;
		
		$nThreads = $1 if /^nThreads=(\d+)$/;
		$nQueries = $1 if /^nQueries=(\d+)$/;
		$nRequests = $1 if /^nRequests=(\d+)$/;
		
		$exUncultured = $1 if /^exUncultured=([01])$/;
		@searchTaxids = split(/,/, $1) if /^searchTaxids=(.+)$/;
		@ignoreTaxids = split(/,/, $1) if /^ignoreTaxids=(.+)$/;
		push @ignoreTaxa, split(/,/, $1) if /^ignoreTaxa=(.+)$/;
		$taxonUCASE = $1 if /^taxonUCASE=([01])$/;
		$mergeDuplicates = $1 if /^mergeDuplicates=([01])$/;
		$ignoreParalogs = $1 if /^ignoreParalogs=([01])$/;
		$ignoreSeqRepeats = $1 if /^ignoreSeqRepeats=([01])$/;
		$ignoreSubspecies = $1 if /^ignoreSubspecies=([01])$/;

		$seqBlast = $1 if /^seqBlast=([01])$/;		
		$taxBlast = $1 if /^taxBlast=([01])$/;
		$alnBlast = $1 if /^alnBlast=([01])$/;
		
		@ranks = split (/,/, $1) if /^ranks=(.+)$/;
	}
	close IN;
}
$maxHits = $nHits unless $maxHits;
my %selfTaxids = ();
if ($selfTax){
	@a = split (/,/, $selfTax);
	foreach (@a){
		@b = split (/:/);
		$selfTaxids{$b[0]} = $b[1];
	}
}
if ($identity and $identity < 1){ $identity *= 100; }
if ($coverage and $coverage < 1){ $coverage *= 100; }
@ignoreTaxa = split (/,/, "unknown,uncultured,unidentified,unclassified,environmental,plasmid,vector,synthetic,phage") unless @ignoreTaxa;


## generate Entrez query text ##

unless ($eqText){
	if ($exUncultured){
		$eqText = "all [filter] NOT(environmental samples[filter] OR metagenomes[orgn])";
	}
	if (@searchTaxids){
		for ($i=0; $i<=$#searchTaxids; $i++){
			next unless ($searchTaxids[$i] =~ /^\d+$/);
			$eqText .= " OR " if $eqText;
			$eqText .= "txid$searchTaxids[$i]\[orgn\]";
		}
	}
	if (@ignoreTaxids){
		$eqText .= " NOT (";
		for ($i=0; $i<=$#ignoreTaxids; $i++){
			next unless $ignoreTaxids[$i] =~ /^\d+$/;
			$eqText .= " OR " if $i;
			$eqText .= "txid$ignoreTaxids[$i]\[orgn\]";
		}
		$eqText .= ")";
	}
}


## read query list ##

print "Reading protein sets...\n";
die "\nError: input directory does not exist.\n" unless (-d "$wkDir/input");
opendir (DIR, "$wkDir/input");
@a = readdir(DIR);
closedir DIR;
foreach (@a){
	next if (/^\./);
	next unless -s "$wkDir/input/$_";
	%h = ('file'=>$_, 'pts'=>[], 'done'=>0);
	if (/^(.+)\.[^\.]+$/){ $ins{$1} = {%h}; }
	else { $ins{$_} = {%h}; }
}
foreach my $set (sort keys %ins){

	# guess input format #
	$s = $ins{$set}{'file'};
	my $intype = "";
	open IN, "<$wkDir/input/$s" or die "\nError: $wkDir/input/$s not readable.\n";
	while(<IN>){
		s/\s+$//; next unless $_; next if /^#/;
		unless ($intype){
			if (length($_) >= 12 and substr($_, 0, 12) eq "LOCUS       "){
				$intype = 'GenBank';
				next;
			}elsif (substr($_, 0, 1) eq ">"){
				$intype = 'FASTA';
				$isInSeq = 1 if (index($_, "|") == -1);
			}else{
				$intype = 'plain list';
			}
		}
		if ($intype eq 'GenBank'){
			push (@{$ins{$set}{'pts'}}, $1) if (/^\s+\/protein_id="(.+)"$/);
		}elsif ($intype eq 'FASTA'){
			if (s/^>//){
				if ($isInSeq){ # user-defined FASTA
					@a = split (/\s+/);
					push (@{$ins{$set}{'pts'}}, $a[0]);
					$inSeqs{$a[0]} = "";
					$t = $a[0];
				}else{ # NCBI FASTA title
					@a = split (/\|/);
					push (@{$ins{$set}{'pts'}}, $a[3]);
				}
			}else{
				$inSeqs{$t} .= $_ if ($isInSeq);
			}
		}elsif ($intype eq 'plain list'){
			next if /^\d/;
			@a = split (/\s+/);
			push @{$ins{$set}{'pts'}}, $a[0];
		}
	}
	close IN;
	s/\.\d+$// for @{$ins{$set}{'pts'}}; # delete version number
	$i = scalar @{$ins{$set}{'pts'}};
	if ($i){ print "  $set: $i protein entries (format: $intype".(", user-defined" x ($isInSeq)).")\n"; }
	else{ print "  $set: no protein entries. Skipped.\n"; }
}
$nPt += scalar @{$ins{$_}{'pts'}} for (keys %ins);
print "Done. $nPt proteins from ".(scalar keys %ins)." sets to BLAST.\n";


## read BLAST reports from previous runs ##

if (-d "$wkDir/blast"){
	print "Reading BLAST results...";
	foreach my $set (keys %ins){
		next unless -d "$wkDir/blast/$set";
		foreach my $pt (@{$ins{$set}{'pts'}}){
			$nDone ++ if (-s "$wkDir/blast/$set/$pt.bla");
		}
		$ins{$set}{'done'} = 1 if ($nDone == scalar(@{$ins{$set}{'pts'}}));
	}
	print " done. $nDone results found, remaining ".($nPt - $nDone)." proteins to BLAST.\n";
}

if ($nPt-$nDone <= 0){
	print "BLAST reports are available for all proteins. blaster.pl exits.\n";
	if ($interactive){
		print "Proceed with HGT prediction (analyer.pl) (yes/NO)?";
		$s = <STDIN>; chomp $s;
		if ($s =~ /^yes$/i or $s =~ /^y$/i){
			my $program = $0;
			$program =~ s/blaster.pl$/analyzer.pl/;
			if (-e $program){
				exec ($^X, $program, $wkDir);
				exit 0;
			}else{
				print "Error: analyzer.pl is not found.\n";
				exit 1;
			}
		}
	}else{ exit 0; }
}


## detect pre-computed BLAST results ##

#################################################################
## This function was inspired by Conor Meehan (cmeehan@itg.be) ##
#################################################################

if ($blastMode and $preBlast and -d $preBlast){
	$i = 0;
	opendir (DIR, $preBlast);
	@a = readdir(DIR);
	closedir DIR;
	foreach (@a){
		next if (/^\./);
		next unless -s "$preBlast/$_";
		$s = $_;
		$s = $1 if /^(.+)\.[^\.]+$/;
		if (exists $ins{$s}){
			$ins{$s}{'pre'} = $_;
			$i ++;
		}
	}
	print "Pre-computed BLAST results are found for $i protein sets.\n";
}


## read local taxonomy database ##

if ($blastMode){
	print "Reading local taxonomy database...";
	die "\nError: local taxonomy database is not found under $taxdump.\n" unless (-e "$taxdump/nodes.dmp" and -e "$taxdump/names.dmp");
	open IN, "<$taxdump/nodes.dmp";
	while (<IN>){
		s/\s+$//;
		@a = split (/\s+\|\s+/);
		%h = ('parent', $a[1], 'rank', $a[2]);
		$taxdumps{$a[0]} = {%h};
	}
	close IN;
	open IN, "<$taxdump/names.dmp";
	while (<IN>){
		s/\s+$//;
		next unless (/scientific name\s*\|$/);
		@a = split (/\s+\|\s+/);
		$taxdumps{$a[0]}{'name'} = $a[1] if (exists $taxdumps{$a[0]});
	}
	close IN;
	print " done. ".scalar(keys %taxdumps)." records read.\n";
}


## read GI to TaxID dictionary ##

if ($blastMode and $gi2taxid and -s $gi2taxid){
	print "Reading GI-to-TaxID dictionary...";
	open IN, "<$gi2taxid" or die "\nError: Invalid dictionary file $gi2taxid.\n";
	while (<IN>){
		s/\s+$//; next unless $_;
		@a = split (/\s+/);
		next unless $#a;
		$gi2taxids{$a[0]} = $a[1];
	}
	close IN;
	print " done. ".scalar(keys %gi2taxids)." records read.\n";
}


## read taxonomic information ##

if (-d "$wkDir/taxonomy"){
	print "Reading taxonomy database...";
	if (-e "$wkDir/taxonomy/taxa.db"){
		open IN, "<$wkDir/taxonomy/taxa.db";
		while (<IN>){
			s/\s+$//; next if /^#/; next unless $_;
			@a = split (/\t/);
			next if exists $dbTaxa{$a[0]};
			my %taxon :shared = ('organism',$a[1],'lineage',$a[2]);
			$i = 3; $taxon{$_} = $a[$i++] for (@ranks);
			$dbTaxa{$a[0]} = \%taxon;
		}
		close IN;
	}
	if (-e "$wkDir/taxonomy/ranks.db"){
		open IN, "<$wkDir/taxonomy/ranks.db";
		while (<IN>){
			s/\s+$//; next if /^#/; next unless $_;
			@a = split (/\t/);
			next if exists $dbRanks{$a[0]};
			$dbRanks{$a[0]} = $a[1];
		}
		close IN;
	}
	print " done. ".(scalar keys %dbTaxa)." taxa and ".(scalar keys %dbRanks)." ranks read.\n";
}else{
	mkdir "$wkDir/taxonomy";
}


## collect 'self' information ##

if (-e "$wkDir/taxonomy/self.info"){
	print "Reading taxonomy of input protein sets...";
	my $nRead = 0;
	open IN, "<$wkDir/taxonomy/self.info";
	while (<IN>){
		s/\s+$//; next if /^#/; next unless $_;
		@a = split (/\t/); next unless @a;
		foreach my $set (keys %ins){
			if ($a[0] eq $set){
				$ins{$set}{'accn'} = $a[1];
				$ins{$set}{'taxid'} = $a[2];
				$ins{$set}{'organism'} = $a[3];
				$nRead ++;
				last;
			}
		}
	}
	close IN;
	print " done. $nRead organisms identified.\n";
}
$s = 1;
foreach my $set (keys %ins){
	unless (exists $ins{$set}{'taxid'}){
		$s = 0;
		last;
	}
}
unless ($s){
	if ($isInSeq){
		print "Defining taxonomy of input protein sets...\n";
		foreach my $set (sort keys %ins){
			next if exists $ins{$set}{'taxid'};
			$ins{$set}{'accn'} = 0; ###################
			$i = 0;
			if (exists $selfTaxids{$set}){
				$i = $selfTaxids{$set};
				$ins{$set}{'taxid'} = $i;
			}elsif ($interactive){
				print "Type the NCBI TaxID of $set (0 if not applicable):";
				$i = <STDIN>; chomp $i;
				$ins{$set}{'taxid'} = $i;
			}
			if ($i){
				@a = get_taxonomy (($i));
				print "  $set is $a[0] ($i).\n";
				$ins{$set}{'organism'} = $a[0];
			}else{
				$ins{$set}{'organism'} = "na";
			}
			open OUT, ">>$wkDir/taxonomy/self.info";
			print OUT $set,"\t",$ins{$set}{'accn'},"\t",$ins{$set}{'taxid'},"\t",$ins{$set}{'organism'},"\n";
			close OUT;
		}
	}else{
		print "Identifying taxonomy of input protein sets...\n";
		my $nIdentified = 0;
		my @accns = (); # proteins whose taxonomic information is to be identified
		foreach my $set (keys %ins){
			next if exists $ins{$set}{'taxid'};
			$ins{$set}{'accn'} = $ins{$set}{'pts'}[0];
			push @accns, $ins{$set}{'accn'};
		}
		if ($blastMode){
			my $query = join (",", @accns);
			my @out = `$blastdbcmd -db $dbBlast -entry $query -outfmt \"%a %g %T %t\"`;
			foreach (@out){
				s/\s+$//; @b = split (/\s+/); $b[0] =~ s/\.\d+$//;
				foreach my $set (keys %ins){
					if ($ins{$set}{'accn'} eq $b[0]){
						last if (exists $ins{$set}{'taxid'} and exists $ins{$set}{'organism'});
						$ins{$set}{'taxid'} = $b[2];
						if (/.*\[([^\[\]]+?)\]$/){ $ins{$set}{'organism'} = $1; }
						else{ $ins{$set}{'organism'} = "na"; }
						print "  ".$ins{$set}{'organism'}." (".$ins{$set}{'taxid'}.")\n";
						open OUT, ">>$wkDir/taxonomy/self.info";
						print OUT $set,"\t",$ins{$set}{'accn'},"\t",$ins{$set}{'taxid'},"\t",$ins{$set}{'organism'},"\n";
						close OUT;
						$nIdentified ++;
						last;
					}
				}
			}
		}else{
			my @ids = ();
			for ($i=0; $i<=$#accns; $i+=20){ # 20 is the maximum number of queries per request
				@a = ();
				for (0..19){ push (@a, $accns[$i+$_]) if ($i+$_<=$#accns); }
				my $iRetry = 0;
				while (1){
					$s = get "$eSearchServer?db=protein&term=".join (",", @a);
					last if (defined $s) and ($s =~ /<Count>\d+<\/Count>/s);
					die "\nFailed to retrieve taxonomic information from NCBI.\n" if ($iRetry >= $nRetry);
					$iRetry ++; sleep $delay; next;
				}
				push (@ids, $1) while ($s =~ s/<Id>(\d+)<\/Id>//s);
				sleep $delay;
			}
			die "\nFailed to identify taxonomy of input protein sets.\n" unless @ids;
			my $iRetry = 0;
			while (1){
				$s = get "$eSummaryServer?db=protein&id=".join (",", @ids);
				last if (defined $s) and ($s =~ /<eSummaryResult>/s);
				die "\nFailed to retrieve taxonomic information from NCBI.\n" if ($iRetry >= $nRetry);
				$iRetry ++; sleep $delay; next;
			}
			@a = (); push (@a, $1) while $s =~ s/<DocSum>(.+?)<\/DocSum>//s;
			foreach my $record (@a){
				$record =~ /<Item Name=\"Caption\" Type=\"String\">(\S+?)<\/Item>/; # Caption is the accession number without version number
				my $set = "";
				foreach (keys %ins){ 
					if ($ins{$_}{'accn'} eq $1){ $set = $_; last; }
				}
				next unless $set;
				$record =~ /<Item Name=\"TaxId\" Type=\"Integer\">(\d+?)<\/Item>/;
				$ins{$set}{'taxid'} = $1;
				$record =~ /<Item Name=\"Title\" Type=\"String\">.*\[([^\[\]]+?)\]<\/Item>/;
				$ins{$set}{'organism'} = $1;
				print "  ".$ins{$set}{'organism'}." (".$ins{$set}{'taxid'}.")\n";
				open OUT, ">>$wkDir/taxonomy/self.info";
				print OUT $set,"\t",$ins{$set}{'accn'},"\t",$ins{$set}{'taxid'},"\t",$ins{$set}{'organism'},"\n";
				close OUT;
				$nIdentified ++;
			}
		}
		@a = ();
		foreach my $set (keys %ins){
			next if exists $dbTaxa{$ins{$set}{'taxid'}};
			push (@a, $ins{$set}{'taxid'});
		}
		get_taxonomy @a;
		print "Done. $nIdentified organisms identified.\n";
	}
}


## perform batch BLAST ##

mkdir "$wkDir/blast" unless -d "$wkDir/blast";
my %pre = (); # pre-computed BLAST results
if ($blastMode or $nRequests == 1){
	foreach my $set (sort keys %ins){
		print "Batch BLAST of $set (".(scalar @{$ins{$set}{'pts'}})." queries) started.\n";
		mkdir "$wkDir/blast/$set" unless -d "$wkDir/blast/$set";
		@a = localtime(time); $s = (++$a[4])."/$a[3]/".($a[5]+=1900)." $a[2]:$a[1]:$a[0]";
		open LOG, ">>$wkDir/blast/$set.log";
		print LOG "Program started at $s.\n";
		print LOG "Number of queries: ". (scalar @{$ins{$set}{'pts'}}) .".\n";
		close LOG;
		if ($blastMode){
		
			# local BLAST #

			%pre = ();
			if (exists $ins{$set}{'pre'}){
			
				# read pre-computed BLAST results #
				
				# In a standard tabular output (-outfmt 6), the text fields are:
				# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
				# The script will only consider seven columns:
				# qseqid (0), sseqid (1), pident (2), qstart (8), qend (9), evalue (10) and bitscore (11)
				# If there is the 12th field, it will be the sequence of the match region.
			
				print "Reading pre-computed BLAST results for $set...";
				$i = 0;
				open IN, "<$preBlast/".$ins{$set}{'pre'};
				while (<IN>){
					s/\s+$//; next unless $_; next if /^#/;
					@a = split (/\t/);
					next if ($evalue and $a[10] ne "*" and $evalue < $a[10]); # % evalue cutoff
					next if ($identity and $a[2] ne "*" and $identity > $a[2]); # % identity cutoff
					%h = ('sseqid'=>$a[1], 'pident'=>$a[2], 'evalue'=>$a[10], 'bitscore'=>$a[11], 'qstart'=>$a[6], 'qend'=>$a[7], 'sseq'=>"");
					$h{'sseq'} = $a[12] if ($#a >= 12 and ($seqBlast or $alnBlast));
					$i ++;
					if (exists $pre{$a[0]}){ push @{$pre{$a[0]}}, {%h}; }
					else{ $pre{$a[0]} = [({%h})]; }
				}				
				close IN;
				print " done. $i hits of ".scalar(keys %pre). " queries read.\n";
			}

			my @queries = ();
			foreach my $query (@{$ins{$set}{'pts'}}){
				next if -s "$wkDir/blast/$set/$query.bla";
				push (@queries, $query);
				next unless $nQueries;
				unless (($#queries+1) % $nQueries){
					local_blast (\@queries, $set);
					@queries = ();
				}
			}
			local_blast (\@queries, $set) if @queries;
		}else{
		
			# single-process http BLAST #
			
			foreach my $query (@{$ins{$set}{'pts'}}){
				next if -s "$wkDir/blast/$set/$query.bla";
				my $return = 0;
				print "  BLASTing $query...";
				my $iRetry = 0;
				while (1){
					$return = http_blast ($query, $set);
					last if $return =~ /\/1$/;
					print " failed.";
					print "\n" and last if ($iRetry >= $nRetry);
					print " Retrying...";
					$iRetry ++; sleep $delay;
				}
				open LOG, ">>$wkDir/blast/$set.log";
				if ($return =~ /\/(\d+)\/1$/){ print LOG "$query\t$1\n"; print " done.\n"; }
				else{ print LOG "$query\tfailed\n"; }
				close LOG;
			}
		}
		@a = localtime(time); $s = (++$a[4])."/$a[3]/".($a[5]+=1900)." $a[2]:$a[1]:$a[0]";
		open LOG, ">>$wkDir/blast/$set.log";
		print LOG "Program ended at $s.\n";
		close LOG;
		unlink "$wkDir/blast.seq";
		print "Batch BLAST of $set (".(scalar @{$ins{$set}{'pts'}})." queries) completed.\n";
	}
}else{

	# multi-process http BLAST #
	
	print "$nRequests http BLAST threads are running in parallel.\n";
	
	my %running = (); # accn -> 1
	my %retry = (); # accn -> times of retry
	my %rSets = (); # protein sets that have been started
	my %failed = ();
	
	while ($nDone + (scalar keys %failed) < $nPt){
		if ($n = $nRequests-scalar(threads->list(threads::running))){
			for (1..$n){
				foreach my $set (keys %ins){
					next if $ins{$set}{'done'};
					unless (exists $rSets{$set}){
						print "Batch BLAST of $set (".(scalar @{$ins{$set}{'pts'}})." queries) started.\n";
						$rSets{$set} = 1;
						mkdir "$wkDir/blast/$set" unless -d "$wkDir/blast/$set";
						@a = localtime(time); $s = (++$a[4])."/$a[3]/".($a[5]+=1900)." $a[2]:$a[1]:$a[0]";
						open LOG, ">>$wkDir/blast/$set.log"; print LOG "Program started at $s.\nNumber of queries: ". (scalar @{$ins{$set}{'pts'}}) .".\n"; close LOG;
					}
					my $started = 0;
					foreach my $query (@{$ins{$set}{'pts'}}){
						next if -s "$wkDir/blast/$set/$query.bla";
						next if exists $running{$query};
						next if exists $failed{"$set|$query"};
						my $thr = threads->create(\&http_blast, $query, $set);
						print "  BLAST of $query started.\n";
						$running{$query} = 1;
						$started = 1;
						select (undef, undef, undef, 0.25); # this is to delay 0.25 seconds
						last;
					}
					last if $started;
					print "Batch BLAST of $set (".(scalar @{$ins{$set}{'pts'}})." queries) completed.\n";
					$ins{$set}{'done'} = 1;
					delete $rSets{$set};
					@a = localtime(time); $s = (++$a[4])."/$a[3]/".($a[5]+=1900)." $a[2]:$a[1]:$a[0]";
					open LOG, ">>$wkDir/blast/$set.log"; print LOG "Program ended at $s.\n"; close LOG;
				}
			}
		}
		if (@a = threads->list(threads::joinable)){
			foreach (@a){
				my ($query, $set, $hits, $return) = split (/\//, $_->join());
				delete $running{$query};
				if ($return){
					$nDone ++;
					print "  BLAST of $query completed. $hits hits retrieved.\n";
					open LOG, ">>$wkDir/blast/$set.log";
					print LOG "$query\t$hits\n";
					close LOG;
				}else{
					$retry{$query} = 0 unless exists $retry{$query};
					$retry{$query} ++;
					if ($retry{$query} >= $nRetry){
						delete $retry{$query};
						$failed{"$set|$query"} = 1;
						print "  BLAST of $query failed.\n";
						open LOG, ">>$wkDir/blast/$set.log"; print LOG "$query\tfailed\n"; close LOG;
					}else{
						print "  BLAST of $query failed. Scheduled to retry.\n";
					}
				}
			}
		}
		sleep $delay;
		if ($retryFailed and %failed and ($nDone + (scalar keys %failed) == $nPt)){
			$retryFailed = 0;
			%failed = ();
			$ins{$_}{'done'} = 0 for (keys %ins);
		}
	}
}

print "Batch BLAST task completed.";


## The following codes have been replaced by a section in HGTector.pl ##

exit 0;
$i = 1;
if ($interactive){
	print "Check completeness of BLAST reports using validator.pl (YES/no)? ";
	$i = 0; $s = <STDIN>; chomp $s;
	$i = 1 if (not $s or $s =~ /^yes$/i or $s =~ /^y$/i);
}
if ($i){
	my $program = $0;
	$program =~ s/blaster.pl$/validator.pl/;
	unless (-e $program){
		print "Error: validator.pl is not found.\n";
		exit 1;
	}
	system ($^X, $program, $wkDir);
}

$i = 1;
if ($interactive){
	print "Restart blaster.pl to fill any missing results (YES/no)? ";
	$i = 0; $s = <STDIN>; chomp $s;
	$i = 1 if (not $s or $s =~ /^yes$/i or $s =~ /^y$/i);
}
exec ($^X, $0, @ARGV) if ($i);
exit 0;


## BLAST using the standalone ncbi-blast+ program ##

sub local_blast {

	# generate query information #
	my ($refQuery, $set) = @_;
	my @queries = @$refQuery;
	$s = "  BLASTing ";
	$s = "  Importing " if ($preBlast and exists $ins{$set}{'pre'});
	print $s.join(",", @queries)."...";
	unlink "$wkDir/blast.seq" if -e "$wkDir/blast.seq";
	my %lengths = ();
	foreach my $query (@queries){
		open (OUT, ">$wkDir/blast/$set/$query.bla");
		print OUT "#NEXUS\nBEGIN QUERY;\n";
		if ($isInSeq){
			$lengths{$query} = length($inSeqs{$query});
			print OUT "\tName=$query;\n\tLength=$lengths{$query};\n\tProduct=na;\n\tOrganism=".$ins{$set}{'organism'}.";\nEND;\n\n";
			open TMP, ">>$wkDir/blast.seq"; print TMP ">".$query."\n".$inSeqs{$query}."\n"; close TMP;
		}else{ # look up query
		
			##### BAD MOVE #####
			my @out = `$blastdbcmd -db $dbBlast -entry $query -outfmt \"%a %g %l %T %s %t\"`;
			foreach (@out){
				s/\s+$//; @b = split (/\s+/);
				last if ($b[0] eq $query or $b[0] =~/^$query\.\d+$/);
			}
			$lengths{$query} = $b[2];
			print OUT "\tGI=$b[1];\n\tAccession=$b[0];\n\tLength=$b[2];\n";
			$s = join (" ", @b[5..$#b]);
			$s =~ /^\s*(.+\S)\s*\[(.+)\]$/; ##### may need error treatment
			print OUT "\tProduct=$1;\n\tOrganism=$2;\nEND;\n\n";
			open TMP, ">>$wkDir/blast.seq"; print TMP ">".$query."\n".$b[4]."\n"; close TMP;
		}
		close OUT;
	}

	my %query2hits = ();
	my %gi2accn = ();
	if ($preBlast and exists $ins{$set}{'pre'}){
	
		# parse pre-computed BLAST results
		foreach my $query (@queries){
			next unless exists $pre{$query};
			@a = @{$pre{$query}};
			for ($i=0; $i<=$#a; $i++){
				next if ($coverage and $coverage > ($a[$i]{'qend'}-$a[$i]{'qstart'}+1)/$lengths{$query}*100); # % coverage cutoff
				@b = split(/;/, $a[$i]{'sseqid'}); # actually, 'sseqid' is only one entry. 'sallseqid' can be multiple.
				my @accns = (); my @gis = ();
				foreach (@b){
					if (/^gi\|(\d+)\|.+\|(.+)\.\d+\|$/){
						$gi2accn{$1} = $2;
						push (@gis, $1);
						push (@accns, $2);
					}else{
						push (@gis, 0);
						push (@accns, $_); # the sseqid is the actual query name (self BLAST)
					}
				}
				# the report contains: subject GIs (all), subject accessions (all), E-value, bit score, percent identity, qstart, qend, aligned part of subject sequence
				$s = join("\t", (join(";", @gis), join(";", @accns), $a[$i]{'evalue'}, $a[$i]{'bitscore'}, $a[$i]{'pident'}, $a[$i]{'qstart'}, $a[$i]{'qend'}, $a[$i]{'sseq'}));
				if (exists $query2hits{$query}){
					push (@{$query2hits{$query}}, $s);
				}else{
					$query2hits{$query} = [ $s ];
				}
			}
		}
	}else{

		# run local BLASTP #
		$s = "$blastp -query $wkDir/blast.seq -db $dbBlast -evalue $evalue -max_target_seqs $nHits -outfmt \"6 qseqid sallgi sallacc evalue bitscore pident qstart qend sseq\"";
		$s .= " -num_threads $nThreads" if ($nThreads > 1);
		my @out = `$s`;
		unlink "$wkDir/blast.seq";
		foreach (@out){
			s/\s+$//; @a = split (/\t/);
			next if ($identity and $identity > $a[5]); # % identity cutoff
			next if ($coverage and $coverage > ($a[7]-$a[6]+1)/$lengths{$a[0]}*100); # % coverage cutoff
			my @gis = split (/;/, $a[1]);
			my @accns = split (/;/, $a[2]);
			s/\.\d+$// for @accns;
			next unless (@gis and @accns);
			foreach ($i=0; $i<=$#gis; $i++){
				$gi2accn{$gis[$i]} = $accns[$i] unless exists $gi2accn{$gis[$i]};
			}
			if (exists $query2hits{$a[0]}){
				push (@{$query2hits{$a[0]}}, join("\t",@a[1..$#a]));
			}else{
				$query2hits{$a[0]} = [ join("\t",@a[1..$#a]) ];
			}
		}
	}

	# get TaxIDs for proteins #
	my %accn2taxid = ();
	if ($gi2taxid){
		if (%gi2taxids){
			foreach my $gi (keys %gi2accn){
				$accn2taxid{$gi2accn{$gi}} = $gi2taxids{$gi} if exists $gi2taxids{$gi};
			}
		}
	}
	if ((scalar keys %accn2taxid) < (scalar keys %gi2accn)){
		open TMP, ">$wkDir/gis.txt";
		foreach my $gi (keys %gi2accn){
			next if exists $accn2taxid{$gi2accn{$gi}};
			print TMP "gi|$gi\n";
		}
		close TMP;
		my @out2 = `$blastdbcmd -dbtype=prot -db $dbBlast -entry_batch $wkDir/gis.txt -target_only -outfmt \"%g %a %T\"`;
		foreach (@out2){
			s/\s+$//; next unless $_;
			@a = split (/\s+/); next if $#a < 2;
			next unless exists $gi2accn{$a[0]};
			next unless $a[2];
			$a[1] =~ s/\.\d+$//;
			$accn2taxid{$a[1]} = $a[2];
		}
		unlink "$wkDir/gis.txt";
	}
	
	# ignore invalid organism names
	my @unwanted = ();
	foreach my $accn (keys %accn2taxid){
		my $taxid = $accn2taxid{$accn};
		unless (exists $taxdumps{$taxid}){ push (@unwanted, $accn); next; }
		my $organism = $taxdumps{$taxid}{'name'};
		if ($organism =~ /^Unresolved/){ push (@unwanted, $accn); next; }
		if ($taxonUCASE and ($organism !~ /^[A-Z]/)){ push (@unwanted, $accn); next; }
		if (@ignoreTaxa){
			foreach (@ignoreTaxa){
				if ($organism =~ /$_/){ push (@unwanted, $accn); last; }
			}
		}
	}
	delete $accn2taxid{$_} for @unwanted;
	foreach (@queries){
		$accn2taxid{$_} = $ins{$set}{'taxid'};
	}

	foreach my $query (@queries){

		# extract valid hits from BLAST result #
		my @hits;
		my %taxids = (); # TaxIDs to look up
		my $isQueryIn = 0; # whether query is among subjects
		if (exists $query2hits{$query}){
			my %usedTaxids = ();
			foreach (@{$query2hits{$query}}){ # each element is a hit line
				@a = split (/\t/);
				# the report contains: subject GIs (all), subject accessions (all), E-value, bit score, percent identity, qstart, qend, aligned part of subject sequence
				my %hit = ('expect'=>$a[2], 'score'=>$a[3], 'identity'=>$a[4], 'coverage'=>sprintf("%.2f", ($a[6]-$a[5]+1)/$lengths{$query}*100), 'sequence'=>$a[7]);
				my @accns = split (/;/, $a[1]); # multiple accns per hit
				s/\.\d+$// for @accns;
				my %taxid2accns = ();
				foreach my $accn (@accns){
					next unless exists $accn2taxid{$accn};
					unless ($isInSeq or $isQueryIn){
						$isQueryIn = 1 if ($accn eq $query);
					}
					if (exists $taxid2accns{$accn2taxid{$accn}}){
						$taxid2accns{$accn2taxid{$accn}} .= "/$accn";
					}else{
						$taxid2accns{$accn2taxid{$accn}} = $accn;
					}
				}
				foreach my $taxid (keys %taxid2accns){
					next if exists $usedTaxids{$taxid};
					$hit{'organism'} = $taxdumps{$taxid}{'name'};
					$hit{'accn'} = order_accns($taxid2accns{$taxid}, $query);
					$hit{'taxid'} = $taxid;
					$taxids{$taxid} = 1;
					$usedTaxids{$taxid} = 1;
					push (@hits, {%hit});
					last if (scalar @hits >= $nHits);
				}
				last if (scalar @hits >= $nHits);
			}
		}

		# perform self BLAST, in case not targetted in previous steps #
		unless ($isQueryIn){
			@a = self_blast $query;
			return "$query/$set/0/0" if (@a == (0,0,0,0));
			my %hit = ('accn'=>$query, 'expect'=>$a[0], 'score'=>$a[1], 'identity'=>$a[2], 'coverage'=>'100.00', 'taxid'=>$ins{$set}{'taxid'}, 'organism'=>$ins{$set}{'organism'}, 'sequence'=>$inSeqs{$query});
			unshift (@hits, {%hit});
			pop @hits if (scalar @hits > $nHits);
		}

		# retrieve complete taxonomy information from local database
		get_taxonomy (keys %taxids) if (%taxids);

		# mark redundant hits from (multiple subspecies or strains of) one species
		if ($ignoreSubspecies){
			my %speciesDB = ();
			for ($i=0; $i<=$#hits; $i++){
				unless (exists $dbTaxa{$hits[$i]{'taxid'}} and exists $dbTaxa{$hits[$i]{'taxid'}}{'species'} and $dbTaxa{$hits[$i]{'taxid'}}{'species'}){
					$hits[$i]{'delete'} = 1; next;
				}
				$s = $dbTaxa{$hits[$i]{'taxid'}}{'species'};
				if (exists $speciesDB{$s}){
					next if ($hits[$i]{'taxid'} == $ins{$set}{'taxid'});
					next if (($hits[$i]{'accn'} eq $query) or ($hits[$i]{'accn'} =~ /^$query\//) or ($hits[$i]{'accn'} =~ /\/$query$/) or ($hits[$i]{'accn'} =~ /\/$query\//));
					$hits[$i]{'delete'} = 1;
				}else{ $speciesDB{$s} = $dbRanks{$s}; }
			}
		}

		# output result
		my ($ntax, $nchar) = (0, 0);
		open (OUT, ">>$wkDir/blast/$set/$query.bla");
		print OUT "BEGIN ORGANISM;\n";
		print OUT "[Accession\tOrganism\tTaxID\tBit-score\tE-value\t\%Identity\t\%Coverage]\n";
		for ($i=0; $i<=$#hits; $i++){
			next if exists $hits[$i]{'delete'};
			$ntax ++;
			$nchar = length($hits[$i]{'sequence'}) if ($hits[$i]{'sequence'} and $nchar < length($hits[$i]{'sequence'}));
			print OUT $hits[$i]{'accn'}."\t".$hits[$i]{'organism'}."\t".$hits[$i]{'taxid'}."\t".$hits[$i]{'score'}."\t".$hits[$i]{'expect'}."\t".$hits[$i]{'identity'}."\t".$hits[$i]{'coverage'}."\n";
			last if $i > $maxHits;
		}
		print OUT ";\nEND;\n\n";
		if ($seqBlast or $alnBlast){
			print OUT "BEGIN DATA;\n";
			print OUT "\tDIMENSIONS NTAX=$ntax NCHAR=$nchar;\n\tFORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\n\tMATRIX\n";
			for ($i=0; $i<=$#hits; $i++){
				next if exists $hits[$i]{'delete'};
				@a = split (/\//, $hits[$i]{'accn'});
				print OUT $a[0]."\t".$hits[$i]{'sequence'}."\n";
				last if $i > $maxHits;
			}
			print OUT ";\nEND;\n\n";
		}
		close OUT;
		open LOG, ">>$wkDir/blast/$set.log";
		print LOG "$query\t$ntax\n";
		close LOG;
	}
	print " done.\n";
	# return "$query/$set/".scalar(@hits)."/1";
}


## BLAST via http connection to NCBI server ##

sub http_blast{ # return: query/set/number of hits/whether successful (0 - failed, 1- successful) # The perl thread join function has some problems, that why I used this inconvenient way.

	# generate query information
	my ($query, $set) = ($_[0], $_[1]);
	my %self = ('name'=>"", 'gi'=>0, 'accession'=>"", 'length'=>0, 'product'=>"", 'organism'=>$ins{$set}{'organism'});
	if ($isInSeq){
		$self{'name'} = $query;
		$self{'length'} = length($inSeqs{$query});
		$self{'product'} = "na";
	}

	# send BLAST request #
	my $isError = 0;
	my $args = "CMD=Put&PROGRAM=blastp&DATABASE=$dbBlast&EXPECT=$evalue&FILTER=m S";
	$args .= "&MAX_NUM_SEQ=$nHits" if ($nHits - 100);
	$args .= "&EQ_TEXT=$eqText" if $eqText;
	my $querykey = $query;
	$querykey = $inSeqs{$query} if $isInSeq; # user-defined query sequence
	$args .= "&QUERY=$querykey";
	# $args = substr($args,0,2048) if (length($args) > 2048); # In some situations, the max size of a URL is 2048 bytes
	my $req = new HTTP::Request POST => $blastServer;
	$req->content_type('application/x-www-form-urlencoded');
	$req->content($args);
	my $response = $ua->request($req);
	my $rid;
	if ($response->content =~ /^    RID = (.*$)/m){ $rid = $1; }
	else{ return "$query/$set/0/0"; } ##########################
	if ($response->content =~ /\s\((\d+) letters\)/){ $self{'length'} = $1; }
	if ($response->content =~ /^    RTOE = (.*$)/m){
		my $second = $1;
		if ($second and $second =~ /^\d+$/){ sleep $second; }else{ sleep $delay; }
	}else{ return "$query/$set/0/0"; }
	while (1){
		sleep $delay;
		$args = "$blastServer?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
		# $args .= "&MAX_NUM_SEQ=$nHits" if ($nHits - 100);
		$req = new HTTP::Request GET => $args;
		$response = $ua->request($req);
		$s = $response->content;
		if ($s =~ /\s+Status=WAITING/m){ next; }
		if ($s =~ /\s+Status=FAILED/m){ $isError = 1; last; }
		if ($s =~ /\s+Status=UNKNOWN/m){ $isError = 1; last; }
		if ($s =~ /\s+Status=READY/m){
			if ($s =~ /\s+ThereAreHits=yes/m){ last; }
			else{ last; } # no hits;
		}
		$isError = 1;
		last;
	}
	if ($isError){ return "$query/$set/0/0"; }

	# retrieve tabular report #
	$args = "$blastServer?CMD=Get&ALIGNMENT_VIEW=Tabular&FORMAT_TYPE=Text&RID=$rid";
	$args .= "&MAX_NUM_SEQ=$nHits&DESCRIPTIONS=$nHits" if ($nHits - 100);
	$req = new HTTP::Request GET => $args;
	$response = $ua->request($req);
	$s = $response->content;
	if (($s !~ /# blastp/) or ($s !~ /# Query:\s/)){ return "$query/$set/0/0"; }
	my %hits = ();
	my $hitid = 0;
	foreach (split(/\n/, $s)){
		if (/^# Query/ and not $isInSeq){
			s/>.*//;
			@a = split (/\|/);
			$self{'gi'} = $a[1];
			$self{'accession'} = $a[3];
			$a[4] =~ /^\s*(.+\S)\s*\[(.+)\]/;
			$self{'product'} = $1;
			$self{'organism'} = $2;
		}
		next if (/^#/);
		@a = split (/\t/);
		next if ($#a < 12);
		next unless $a[12];
		next if ($identity and $identity > $a[2]); # % identity cutoff
		@b = split (/;/, $a[1]);
		$hitid ++;
		foreach (@b){
			next unless (/\|$/);
			@c = split (/\|/);
			$c[3] =~ s/\.\d+$//;
			my %hit = ('id', $hitid, 'accn', $c[3], 'expect', $a[11], 'score', $a[12], 'identity', $a[2], 'qlength', $a[8]-$a[7]+1);
			$hits{$c[1]} = {%hit};
		}
	}

	# retrieve taxonomy information #
	$i = 0; # count
	@a = (); # all results
	@b = (); # subset of GIs
	foreach (keys %hits){
		$i ++;
		push (@b, $_);
		if ($i == 190){ # in some situations, a URI should not exceed ~2000 characters, which is approximately 190 GIs.
			while (1){
				$s = get $eSummaryServer."?db=protein&id=".join (",", @b);
				last if (defined $s);
				sleep $delay;
			}
			push (@a, $1) while ($s =~ s/<DocSum>(.+?)<\/DocSum>//s);
			$i = 0; @b = ();
		}
	}
	if (@b){
		while (1){
			$s = get $eSummaryServer."?db=protein&id=".join (",", @b);
			last if (defined $s);
			sleep $delay;
		}
		push (@a, $1) while ($s =~ s/<DocSum>(.+?)<\/DocSum>//s);
	}

	foreach (@a){
		/<Id>(\d+)<\/Id>/;
		$s = $1;
		/<Item Name=\"TaxId\" Type=\"Integer\">(\d+?)<\/Item>/;
		$hits{$s}{'taxid'} = $1;
		/<Item Name=\"Title\" Type=\"String\">.*\[([^\[\]]+?)\]<\/Item>/;
		$hits{$s}{'organism'} = $1;
		/<Item Name=\"Length\" Type=\"Integer\">(\d+)<\/Item>/;
		$hits{$s}{'length'} = $1;
	}

	# discard hits whose taxonomy information is unidentified #
	foreach (keys %hits){
		delete $hits{$_} and next unless exists $hits{$_}{'taxid'};
		delete $hits{$_} unless $hits{$_}{'taxid'};
		delete $hits{$_} and next unless exists $hits{$_}{'organism'};
		delete $hits{$_} unless $hits{$_}{'organism'};
		delete $hits{$_} if ($hits{$_}{'organism'} =~ /^Unresolved/);
		delete $hits{$_} if ($taxonUCASE and ($hits{$_}{'organism'} !~ /^[A-Z]/));
	}

	# discard hits whose organism names contain specified strings #
	if (@ignoreTaxa){
		foreach (keys %hits){
			foreach $s (@ignoreTaxa){
				if ($hits{$_}{'organism'} =~ /$s/){
					delete $hits{$_};
					last;
				}
			}
		}
	}

	# merge identical proteins (defined by NCBI) from one organism #
	$i = 0; # current hit ID
	$j = 0; # current taxid
	$s = ""; # current key
	foreach ( sort {$hits{$a}{'id'} <=> $hits{$b}{'id'}} keys %hits){
		if ($hits{$_}{'id'} != $i){
			$i = $hits{$_}{'id'};
			$j = $hits{$_}{'taxid'};
			$s = $_;
		}else{
			if ($hits{$_}{'taxid'} == $j){
				$hits{$s}{'accn'} .= "/".$hits{$_}{'accn'};
				delete $hits{$_};
			}
		}
	}

	# merge duplicated hits (proteins with same bit score from one organism) #
	if ($mergeDuplicates){
		$i = 0; # current bit score
		$j = 0; # current taxid
		$s = ""; # current key
		foreach ( sort {($hits{$b}{'score'} <=> $hits{$a}{'score'}) or ($hits{$a}{'taxid'} <=> $hits{$b}{'taxid'})} keys %hits){
			if (($hits{$_}{'score'} != $i) or ($hits{$_}{'taxid'} != $j)){
				$i = $hits{$_}{'score'};
				$j = $hits{$_}{'taxid'};
				$s = $_;
			}else{
				$hits{$s}{'accn'} .= "/".$hits{$_}{'accn'};
				delete $hits{$_};
			}
		}
	}

	# reorder accession numbers #
	foreach (keys %hits){
		$hits{$_}{'accn'} = order_accns ($hits{$_}{'accn'}, $query);
	}

	# retrieve taxonomy information from NCBI server #
	$i = 0; @a = ();
	foreach (keys %hits){
		next if exists $dbTaxa{$hits{$_}{'taxid'}};
		push @a, $hits{$_}{'taxid'};
		$i ++;
		if ($i == 150){
			get_taxonomy @a;
			$i = 0; @a = (); sleep $delay;
		}
	}
	get_taxonomy @a if @a;

	# check whether BLAST result contains query itself #
	my $isQueryIn = 0;
	unless ($isInSeq){
		foreach (keys %hits){
			@a = split (/\//, $hits{$_}{'accn'});
			foreach $s (@a){ if ($query eq $s){ $isQueryIn = $hits{$_}{'id'}; last; }}
			if ($isQueryIn){ $self{'length'} = $hits{$_}{'length'} unless $self{'length'}; last; }
		}
	}

	# find out length of query sequence, in case not in previous steps #

	unless ($self{'length'} or $isInSeq or $blastMode){
		while (1){
			$s = get "$eSearchServer?db=protein&term=$query";
			last if (defined $s);
			sleep $delay;
		}
		$s =~ /<Id>(\d+)<\/Id>/;
		while (1){
			$s = get "$eSummaryServer?db=protein&id=$1";
			last if (defined $s);
			sleep $delay;
		}
		$s =~ /<Item Name=\"Length\" Type=\"Integer\">(\d+)<\/Item>/;
		$self{'length'} = $1;
	}
	next unless $self{'length'};
	$hits{$_}{'coverage'} = sprintf("%.2f", $hits{$_}{'qlength'}/$self{'length'}*100) for (keys %hits);

	# perform self BLAST, in case not targetted in previous steps #

	unless ($isQueryIn){
		@a = self_blast $query;
		return "$query/$set/0/0" if (@a == (0,0,0,0));
		my %hit = ('id', 0, 'accn', $query, 'expect', $a[0], 'score', $a[1], 'identity', $a[2], 'coverage', '100.00', 'taxid', $ins{$set}{'taxid'}, 'organism', $ins{$set}{'organism'});
		if ($self{'length'}){ $hit{'length'} = $self{'length'}; }
		elsif ($a[3]){ $hit{'length'} = $a[3]; }
		else{ $hit{'length'} = 0; }
		$hits{$query} = {%hit};
	}

	# in multiple hits from one species, only keep the hit with the highest bit score #
	if ($ignoreSubspecies){
		my %speciesDB = ();
		foreach (sort {($hits{$b}{'score'} <=> $hits{$a}{'score'}) or ($hits{$a}{'id'} <=> $hits{$b}{'id'})} keys %hits){
			unless (exists $hits{$_}{'taxid'} and exists $dbTaxa{$hits{$_}{'taxid'}} and exists $dbTaxa{$hits{$_}{'taxid'}}{'species'} and $dbTaxa{$hits{$_}{'taxid'}}{'species'}){
				delete $hits{$_}; next;
			}
			$s = $dbTaxa{$hits{$_}{'taxid'}}{'species'};
			unless (exists $speciesDB{$s}){
				$speciesDB{$s} = $dbRanks{$s};
			}else{
				next if ($hits{$_}{'taxid'} == $ins{$set}{'taxid'});
				next if (($hits{$_}{'accn'} eq $query) or ($hits{$_}{'accn'} =~ /^$query\//) or ($hits{$_}{'accn'} =~ /\/$query$/) or ($hits{$_}{'accn'} =~ /\/$query\//));
				delete $hits{$_};
			}
		}
	}

	# remove excessive hits #

	if ($maxHits){
		$i = 0;
		foreach (sort {$hits{$b}{'score'} <=> $hits{$a}{'score'}} keys %hits){
			$i ++;
			delete $hits{$_} if ($i > $maxHits);
		}
	}

	# create output file #

	open (OUT, ">$wkDir/blast/$set/$query.bla");
	print OUT "#NEXUS\nBEGIN QUERY;\n";
	if ($isInSeq){ print OUT "\tName=".$self{'name'}.";\n\tLength=".$self{'length'}.";\n\tProduct=na;\n\tOrganism=".$ins{$set}{'organism'}.";\nEND;\n\n"; }
	else{ print OUT "\tGI=".$self{'gi'}.";\n\tAccession=".$self{'accession'}.";\n\tLength=".$self{'length'}.";\n\tProduct=".$self{'product'}.";\n\tOrganism=".$self{'organism'}.";\nEND;\n\n"; }

	# retrieve taxonomy report (using TaxBLAST)

	if ($taxBlast){
		$args = "$blastServer?CMD=Get&FORMAT_TYPE=HTML&FORMAT_OBJECT=TaxBlast&RID=$rid&ALIGNMENTS=100000";
		$req = new HTTP::Request GET => $args;
		$response = $ua->request($req);
		$s = $response->content;
		if (($s !~ /Tax BLAST Report/) or ($s !~ /Lineage Report/)){
			print OUT "BEGIN ERROR;\nRetrieval of taxonomy report failed.\n;\nEND;\n\n";
		}else{
			$t = 0; # reading status
			foreach (split(/\n/, $s)){
				if ($_ eq "<B><A NAME=lineage>Lineage Report</A></B><BR>"){
					print OUT "BEGIN LINEAGE;\n";
					$t = 1; next;
				}
				if (($_ eq "</FONT></PRE><HR>") and $t){
					print OUT ";\nEND;\n\n";
					$t = 0; next;
				}
				if ($t){
					s/\<.*?>//g;
					s/( hit[ s] \[.*?\])(.*)$/$1/;
					if (@ignoreTaxa){
						$i = 0;
						foreach $s (@ignoreTaxa){
							if (/$s/){ $i = 1; last; }
						}
						next if $i;
					}
					print OUT "$_\n";
					next;
				}
			}
		}
	}

	# output hit table #

	print OUT "BEGIN ORGANISM;\n";
	print OUT "[Accession\tOrganism\tTaxID\tBit-score\tE-value\t\%Identity\t\%Coverage]\n";
	foreach (sort {($hits{$b}{'score'} <=> $hits{$a}{'score'}) or ($hits{$a}{'id'} <=> $hits{$b}{'id'})} keys %hits){
		print OUT $hits{$_}{'accn'}."\t".$hits{$_}{'organism'}."\t".$hits{$_}{'taxid'}."\t".$hits{$_}{'score'}."\t".$hits{$_}{'expect'}."\t".$hits{$_}{'identity'}."\t".$hits{$_}{'coverage'}."\n";
	}
	print OUT ";\nEND;\n\n";

	# retrieve hit sequences #

	if ($seqBlast){
		$i = 0; # count
		$s = ""; # all results
		@b = (); # subset of GIs
		foreach (keys %hits){
			$i ++;
			push (@b, $_);
			if ($i == 190){ # a URI should not exceed ~2000 characters, which is approximately 190 GIs.
				while (1){
					$t = get $eFetchServer."?db=protein&rettype=FASTA&id=".join (",", @b);
					last if (defined $t);
					sleep $delay;
				}
				$s .= $t;
				$i = 0; @b = ();
			}
		}
		if (@b){
			while (1){
				$t = get $eFetchServer."?db=protein&rettype=FASTA&id=".join (",", @b);
				last if (defined $t);
				sleep $delay;
			}
			$s .= $t;
		}
		$i = 0; # current GI
		foreach (split (/\n/, $s)){
			next unless $_;
			if (/^>gi\|(\d+)\|/){
				$i = $1;
				$hits{$i}{'sequence'} = "";
			}else{
				$hits{$i}{'sequence'} .= $_ if $i;
			}
		}
	}

	# retrieve multiple sequence alignment (conflicts seqBlast)

	if ($alnBlast and not $seqBlast){
		$args = "$blastServer?CMD=Get&ALIGNMENT_VIEW=FlatQueryAnchoredNoIdentities&FORMAT_TYPE=Text&RID=$rid&ALIGNMENTS=100000";
		$req = new HTTP::Request GET => $args;
		$response = $ua->request($req);
		$s = $response->content;
		if (($s !~ /blastp/i) or ($s !~ /\nQuery=\s/)){
			print OUT "BEGIN ERROR;\nRetrieval of multiple sequence alignment failed.\n;\nEND;\n\n";
		}else{
			@c = split(/\n/, $s);
			$t = 0; # reading status
			my $iBlock = -1; # block ID
			my %seqs; # sequence alignment
			foreach (@c){
				if ($_ eq "ALIGNMENTS"){ $t = 1; next; }
				if ($t and /^\s/){ $t = 0; last; }
				next unless $t;
				unless ($_){ $iBlock ++; next; } # Start a new block if empty line
				@a = split(/\s+/,substr($_,0,19)); # id
				#next if ($a[0] eq "Query");
				$_ =~ /(\s\s\d*$)/;
				$s = substr($_,19,length($_)-19-length($1)); # sequence
				$s =~ s/\s/-/g;
				$seqs{$a[0]} = "-"x($iBlock*60) unless (exists $seqs{$a[0]}); # add new sequence
				$seqs{$a[0]} .= $s; # add new sequence
			}
			$i = 0;
			foreach $s(keys %seqs){
				$i = length($seqs{$s}) if (length($seqs{$s}) > $i);
			}
			foreach $s(keys %seqs){
				$seqs{$s} .= "-"x($i - length($seqs{$s})) if (length($seqs{$s}) < $i);
			}
			foreach $t (keys %hits){
				@a = split /\//, $hits{$t}{'accn'};
				foreach (@a){
					if (exists $seqs{$_}){
						$hits{$t}{'sequence'} = $seqs{$_};
						last;
					}
				}
			}
		}
	}

	# output sequences #

	if ($seqBlast or $alnBlast){
		$i = 0; # number of sequences
		$j = 0; # maximum length of sequence
		foreach (keys %hits){
			next unless exists $hits{$_}{'sequence'};
			$i ++;
			$j = length ($hits{$_}{'sequence'}) if (length ($hits{$_}{'sequence'}) > $j);
		}
		print OUT "BEGIN DATA;\n";
		print OUT "\tDIMENSIONS NTAX=$i NCHAR=$j;\n\tFORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\n\tMATRIX\n";
		foreach (sort {($hits{$a}{'id'} <=> $hits{$b}{'id'})} keys %hits){
			next unless exists $hits{$_}{'sequence'};
			@a = split (/\//, $hits{$_}{'accn'});
			print OUT $a[0]."\t".$hits{$_}{'sequence'}."\n";
		}		
		print OUT ";\nEND;\n\n";
	}
	close OUT;
	return "$query/$set/".scalar(keys %hits)."/1";
}


## blast the query sequence against itself ##

sub self_blast {

	# parameter: query
	# return: (e-value, bit-score, identity, length)
	# return (0, 0, 0, 0) if failed
	
	my ($query, $length) = ($_[0], 0);
	$length = length($inSeqs{$query}) if ($isInSeq);
	if ($blastMode){ # local mode
		if ($isInSeq){
			open TMP, ">$wkDir/blast.seq"; print TMP ">$query\n".$inSeqs{$query}."\n"; close TMP;
		}else{
			##### BAD MOVE #####
			my @out = `blastdbcmd -db $dbBlast -entry $query -outfmt \"%f\"`;
			my $queryseq = ""; my $isEnd = 0;
			foreach (@out){
				s/\s+$//;
				if (/^>/){ if ($isEnd){ last; }else{ $isEnd = 1; next; } }
				$queryseq .= $_;
			}
			$length = length($queryseq);
			open TMP, ">$wkDir/blast.seq"; print TMP ">$query\n$queryseq\n"; close TMP;
		}
		my @out = `$blastp -query $wkDir/blast.seq -subject $wkDir/blast.seq -outfmt \"6 evalue bitscore pident\"`;
		@a = split (/\t/, $out[0]);
		return (0,0,0,0) if ($#a < 2);
		$a[2] =~ s/\s+$//;
		return ($a[0], $a[1], $a[2], $length);
			
	}else{ # http mode
		my $isError = 0;
		if ($isInSeq){ $s = get "$blastAlignServer?CMD=Put&PROGRAM=blastp&DATABASE=$dbBlast&QUERY=".$inSeqs{$query}."&SUBJECTS=".$inSeqs{$query}; }
		else{ $s = get "$blastAlignServer?CMD=Put&PROGRAM=blastp&DATABASE=$dbBlast&QUERY=$query&SUBJECTS=$query"; };
		return (0,0,0,0) unless (defined $s and $s =~ /^    RID = (.*$)/m);
		my $rid = $1;
		if ($s =~ /\s\((\d+) letters\)/ and not $length){ $length = $1; }
		while (1){
			$s = get "$blastServer?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
			return (0,0,0,0) unless defined $s and $s =~ /\s+Status=(.+)/m;
			if ($1 eq "WAITING"){ sleep $delay; next; }
			if ($1 eq "FAILED" or $1 eq "UNKNOWN"){ $isError = 1; last; }
			if ($1 eq "READY" and $s =~ /\s+ThereAreHits=yes/m){ last; }
			$isError = 1; last;
		}
		return (0,0,0,0) if $isError;
		$s = get "$blastServer?CMD=Get&ALIGNMENT_VIEW=Tabular&FORMAT_TYPE=Text&RID=$rid";
		return (0,0,0,0) unless (defined $s and $s =~ /# blastp/ and $s =~ /# Query:\s/);
		$s =~ /\<PRE\>(.*)\<\/PRE\>/s; $s = $1;
		foreach (split(/\n/, $s)){
			next unless $_;
			next if (/^#/);
			my @a = split (/\t/);
			$a[12] =~ s/^\s+//;
			return ($a[11], $a[12], $a[2], $length);
		}
	}
}


## retrieve taxonomy information from remote or local database ##

sub get_taxonomy{
	
	# parameter: array of TaxIDs
	# return: array of organism names in same order
	# other: write global variable %dbTaxa and %dbRanks, write file taxa.db and ranks.db

	my @organisms = ();
	my %taxa2w = (); # taxa to write
	my %ranks2w = (); # ranks to write
	my %hRanks = map { $_ => 1 } @ranks;

	if ($blastMode){ # local taxonomy database
		foreach my $taxid (@_){
			next if exists $dbTaxa{$taxid};
			if (exists $taxdumps{$taxid}){
				my $id = $taxid;
				my $pid = $taxdumps{$taxid}{'parent'};
				my $name = $taxdumps{$taxid}{'name'};
				my $rank = $taxdumps{$taxid}{'rank'};
				push (@organisms, $name);
				my %taxon :shared = ('organism', $name);
				$taxon{$_} = "" for (@ranks);
				$taxon{$rank} = $id if ($rank and exists $hRanks{$rank});
				while (1){
					last unless ($pid and $pid != $id);
					$name = $taxdumps{$pid}{'name'};
					$rank = $taxdumps{$pid}{'rank'};
					unless (exists $dbRanks{$pid}){
						$dbRanks{$pid} = $name;
						$ranks2w{$pid} = $name;
					}
					$taxon{'lineage'} .= "/$pid";
					$taxon{$rank} = $pid if (exists $hRanks{$rank});
					$id = $pid;
					$pid = $taxdumps{$id}{'parent'};
				}
				$dbTaxa{$taxid} = \%taxon;
				$taxa2w{$taxid} = {%taxon};
			}else{
				unless (exists $badTaxids{$taxid}){
					print " Warning: TaxID $taxid is not found in the database.";
					$badTaxids{$taxid} = 1;
					open OUT, ">>$wkDir/taxonomy/invalid.taxids";
					print OUT $taxid."\n";
					close OUT;
				}
			}
		}

	}else{ # remote taxonomy database
		my $iRetry = 0;
		my @taxids2get = ();
		foreach (@_){ push (@taxids2get, $_) unless exists $dbTaxa{$_}; }
		while (1){
			$s = get "$eFetchServer?db=taxonomy&id=".join (",", @taxids2get);
			last if (defined $s);
			exit "\nFailed to retrieve taxonomic information from NCBI.\n" if ($iRetry >= $nRetry);
			$iRetry ++; sleep $delay; next;
		}
		$s =~ s/<TaxaSet>//;
		while ($s =~ s/\n<Taxon>\s+<TaxId>(\d+)<\/TaxId>\s+<ScientificName>(.+?)<\/ScientificName>(.+?)\n<\/Taxon>//s){
			my $id = $1; $t = $3;
			push (@organisms, $2);
			my %taxon :shared = ('organism', $2);
			$taxon{$_} = "" for (@ranks);
			while ($t =~ s/<Taxon>\s+<TaxId>(\d+)<\/TaxId>\s+<ScientificName>(.+?)<\/ScientificName>\s+<Rank>(.+?)<\/Rank>\s+<\/Taxon>//s){
				unless (exists $dbRanks{$1}){
					$dbRanks{$1} = $2;
					$ranks2w{$1} = $2;
				}
				$taxon{'lineage'} .= "/$1";
				$taxon{$3} = $1 if (exists $hRanks{$3});
			}
			if ($t =~ /<\/ParentTaxId>\s+<Rank>(\S+)<\/Rank>/){ # the taxon itself is a rank
				foreach (@ranks){
					next if $taxon{$_};
					if ($1 eq $_){ $taxon{$_} = $id; last; }
				}
			}
			$dbTaxa{$id} = \%taxon;
			$taxa2w{$id} = {%taxon};
		}
	}
	if (%taxa2w){
		open OUT, ">>$wkDir/taxonomy/taxa.db";
		foreach my $id (keys %taxa2w){
			print OUT $id."\t".$taxa2w{$id}{'organism'}."\t".$taxa2w{$id}{'lineage'};
			print OUT "\t".$taxa2w{$id}{$_} for (@ranks);
			print OUT "\n";
		}
		close OUT;
	}
	if (%ranks2w){
		open OUT, ">>$wkDir/taxonomy/ranks.db";
		foreach my $id (keys %ranks2w){
			print OUT $id."\t".$ranks2w{$id}."\n";
		}
		close OUT;
	}
	return @organisms;
}


# reorder accession numbers #
# order: NP_ > XP_ = YP_ = ZP_ = AP_ > everything else

sub order_accns{
	return $_[0] unless $_[0] =~ /\//;
	my @accns0 = split (/\//, $_[0]);
	my @accns1 =  ();
	foreach (@accns0){
		if (/_/){ unshift (@accns1, $_); }
		else{ push (@accns1, $_); }
	}
	@accns0 = ();
	foreach (@accns1){
		if (/^NP_/){ unshift (@accns0, $_); }
		else{ push (@accns0, $_); }
	}
	@accns1 = ();
	foreach (@accns0){
		if ($#_ and /^$_[1]$/){ unshift (@accns1, $_); }
		else{ push (@accns1, $_); }
	}
	return join ("/", @accns1);
}


