#!/usr/bin/perl

use strict;
use warnings;
use threads;
use threads::shared;
use LWP::Simple;
$| = 1;

print "
-> Searcher: Batch sequence homology searching and filtering. <-
";

print "

Usage:
  Place input files in input/, one file per protein set, then execute:
  perl searcher.pl <working directory>

Input file format:
  1. Multi-Fasta format
  2. List of protein names (must be included in the database)
  Note: GenBank-style names will be parsed.

Output:
  search/sample_name/protein_name.txt
  search/sample_name.log

" and exit unless @ARGV;


## streamline ##
# 1. read and validate configurations
# 2. read input sets and proteins
# 3. read and validate previous search results
# 4. detect pre-computed search results
# 5. read local databases
# 6. read previous taxonomy records
# 7. identify 'self' taxonomy
# 8. run batch homology search and/or import results


## public variables ##

my $s; my @a; my @b; my @c; my %h;


## global variables ##

my $wkDir = $ARGV[0];					# working directory

my %ins = ();							# the master data structure of the whole procedure
										# set => (
										#	file (file)
										#	taxid (str)
										#	organism (str)
										#	prefile (file)
										#	done (boolean)
										#	prots (array...)
										# prot => (
										#	name (accn or user-defined)
 										#	gi, accn (only for GenBank-style sequence titles)
										#	product (after the 1st whitespace)
										#	seq (necessary for local searches, optional for remote BLAST)
										#	done (boolean)
										#	hits (array...)
										# hit => (
										#	sseqid, pident, evalue, bitscore, qstart, qend
										#	sseq (only for BLAST)
										#	name (accn or full name)

my $nProt = 0;							# total number of proteins
my $nDone = 0;							# number of proteins that have been searched already

my $retryFailed = 1;

my %inSeqs = ();						# query protein sequences (if applicable)
my %taxdumps = ();						# local taxonomy database: taxid => (name, parent, rank)
my %prot2taxids = ();					# protein to TaxID dictionary (the protein name may be GI or accn)
my %dbTaxa :shared = ();				# taxa.db
my %dbRanks :shared = ();				# ranks.db
my %badTaxids = ();						# TaxIDs that don't exist in the local taxonomy database


## subroutines ##

sub http_blast;							# parameters: query (accn. no.), set; return: 1 - succeeded, 0 - fail
sub local_search;						# same as above
sub self_align;							# search against itself
sub get_taxonomy;						# paramter: array of TaxIDs
sub stem_name;							# get stem file name
sub seq_title;							# parse sequence title
sub order_accns;						# reorder accession number


## program parameters ##

my $interactive = 1;					# interactive or automatic mode
my $searchTool = "blast";				# protein sequence similarity search tool
my $preSearch = "";						# directory of pre-computed search results.
my $selfTax = "";						# taxonomy of input protein sets

# databases
my $protdb = "";						# protein database for homolgy search
my $taxdump = "";						# directory of the NCBI taxonomy database (nodes.dmp and names.dmp)
my $prot2taxid = "";					# protein name / GI to TaxID dictionary file

# search cutoffs
my $nHits = 500;						# number of hits to return
my $maxHits = 0;						# maximum number of valid hits to preserve. if 0 then = nHits
my $evalue = 1e-5;						# maximum E-value cutoff
my $identity = 0;						# minimum percent identity cutoff
my $coverage = 0;						# minimum query coverage cutoff

# taxonomy filters
my $mergeDuplicates = 1;				# ignore hits with same taxon names and bit scores
my $taxonUCASE = 0;						# ignore taxon names that do not start with a capital letter
my @ignoreTaxa = split (/,/, "unknown,uncultured,unidentified,unclassified,environmental,plasmid,vector,synthetic,phage");							
										# ignore taxon names containing these words
my $ignoreParalogs = 1;					# ignore potential paralogs (hits with same taxon names but different bit scores)
my $ignoreSeqRepeats = 1;				# ignore repeated sequences (hits targetting different regions of the same protein)
my $ignoreSubspecies = 0;				# ignore more than one subspecies from the same species
my @ranks = ('species', 'genus', 'family', 'order', 'class', 'phylum');
										# record the TaxIDs on these ranks for each hit

# search tool behavior
my $threads = 0;						# multiple threads (0 for all CPU cores)
my $queries = 0;						# multiple queries per run (0 for all sequences per sample)
my $getAln = 0;							# retrieve aligned part of subject sequence (for BLAST only)
my $blastdbcmd = "blastdbcmd";
my $blastp = "blastp";
my $rapsearch = "rapsearch";
my $prerapsearch = "prerapsearch";
my $diamond = "diamond";

# remote BLAST behavior
my $httpBlast = 0;						# call the NCBI server to perform BLAST and to retrieve sequence / taxonomy information
my $httpDb = "nr";						# NCBI BLAST database
my $requests = 1;						# number of requests for http BLAST (default: 1)
my $retries = 5;						# maximum number of retries
my $delay = 30;							# time (seconds) between two http requests
my $timeout = 600;						# time (seconds) to give up waiting
my $exUncultured = 1;					# exclude uncultured and environmental samples
my @searchTaxids = ();					# search under the following taxon groups (taxids)
my @ignoreTaxids = ();					# ignore organisms under the following taxids
my $eqText = "";						# entrez query parameter
my $taxBlast = 0;						# retrieve taxonomy report
my $seqBlast = 0;						# retrieve hit sequences
my $alnBlast = 0;						# retrieve multiple sequence alignment (conflicts seqBlast)
my $blastServer = "http://blast.ncbi.nlm.nih.gov/Blast.cgi";
my $blastAlignServer = "http://blast.ncbi.nlm.nih.gov/BlastAlign.cgi";
my $eSearchServer = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
my $eFetchServer = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
my $eSummaryServer = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi";


## read configuration ##

if (-e "$wkDir/config.txt"){
	open IN, "<$wkDir/config.txt";
	while (<IN>){
		s/#.*$//; s/\s+$//g; s/^\s+//g; next unless $_;
		$interactive = $1 if /^interactive=([01])$/;
		$selfTax = $1 if /^selfTax=(.+)$/;
		@ranks = split (/,/, $1) if /^ranks=(.+)$/;
		$preSearch = $1 if /^preSearch=(.+)$/;
		$preSearch =~ s/\/$//;
		
		$searchTool = $1 if /^searchTool=(.+)$/;
		$threads = $1 if /^threads=(\d+)$/;
		$queries = $1 if /^queries=(\d+)$/;
		$requests = $1 if /^requests=(\d+)$/;
		$getAln = $1 if /^getAln=([01])$/;

		$blastdbcmd = $1 if /^blastdbcmd=(.+)$/;
		$blastp = $1 if /^blastp=(.+)$/;
		$rapsearch = $1 if /^rapsearch=(.+)$/;
		$prerapsearch = $1 if /^prerapsearch=(.+)$/;
		$diamond = $1 if /^diamond=(.+)$/;

		$protdb = $1 if /^protdb=(.+)$/;
		$taxdump = $1 if /^taxdump=(.+)$/;
		$taxdump =~ s/\/$//;
		$prot2taxid = $1 if /^prot2taxid=(.+)$/;
		
		$evalue = $1 if /^evalue=(.+)$/;
		$identity = $1 if /^identity=(.+)$/;
		$coverage = $1 if /^coverage=(.+)$/;
		$nHits = $1 if /^nHits=(\d+)$/;
		$maxHits = $1 if /^maxHits=(\d+)$/;		

		$taxonUCASE = $1 if /^taxonUCASE=([01])$/;
		$mergeDuplicates = $1 if /^mergeDuplicates=([01])$/;
		$ignoreParalogs = $1 if /^ignoreParalogs=([01])$/;
		$ignoreSeqRepeats = $1 if /^ignoreSeqRepeats=([01])$/;
		$ignoreSubspecies = $1 if /^ignoreSubspecies=([01])$/;
		@ignoreTaxa = split(/,/, $1) if /^ignoreTaxa=(.+)$/;
		@ignoreTaxa = () if /^ignoreTaxa=$/;
		
		$httpBlast = $1 if /^httpBlast=(\d)$/;
		$httpDb = $1 if /^httpDb=(.+)$/;
		$retries = $1 if /^retries=(\d+)$/;
		$delay = $1 if /^delay=(\d+)$/;
		$eqText = $1 if /^eqText=(.+)$/;
		$exUncultured = $1 if /^exUncultured=([01])$/;
		@searchTaxids = split(/,/, $1) if /^searchTaxids=(.+)$/;
		@ignoreTaxids = split(/,/, $1) if /^ignoreTaxids=(.+)$/;
		$seqBlast = $1 if /^seqBlast=([01])$/;		
		$taxBlast = $1 if /^taxBlast=([01])$/;
		$alnBlast = $1 if /^alnBlast=([01])$/;
		
		$blastServer = $1 if /^blastServer=(.+)$/;
		$blastAlignServer = $1 if /^blastAlignServer=(.+)$/;
		$eSearchServer = $1 if /^eSearchServer=(.+)$/;
		$eFetchServer = $1 if /^eFetchServer=(.+)$/;
		$eSummaryServer = $1 if /^eSummaryServer=(.+)$/;
	}
	close IN;
}


## validate parameters ##

$maxHits = $nHits unless $maxHits;
$identity *= 100 if ($identity and $identity < 1);
$coverage *= 100 if ($coverage and $coverage < 1);

if (lc($searchTool) eq "blast"){ $searchTool = "BLAST"; }
elsif (lc($searchTool) eq "rapsearch"){ $searchTool = "RAPSearch2"; }
elsif (lc($searchTool) eq "rapsearch2"){ $searchTool = "RAPSearch2"; }
elsif (lc($searchTool) eq "diamond"){ $searchTool = "DIAMOND"; }
elsif (lc($searchTool) eq "customized"){ $searchTool = "customized"; }
else{ die "Error: Invalid search tool: $searchTool.\n"; }

if ($searchTool eq "customized" and not $preSearch){ die "Error: You must provide pre-computed search results for \"customized\" search tool.\n"; }
if ($preSearch and not -d $preSearch){ die "Error: Invalid directory for pre-computed search results: $preSearch\n"; }

if ($searchTool ne "BLAST"){ $httpBlast = 0; }
if ($httpBlast and not $protdb){ $protdb = "nr"; }

unless ($threads){ # attempt to get number of CPU cores
	if (-e "/proc/cpuinfo"){
		$threads = `grep -c ^processor /proc/cpuinfo`;
	}elsif ($^O eq "darwin"){
		$threads = `sysctl -n hw.ncpu`;
	}elsif ($^O eq "MSWin32"){
		$threads = `echo %NUMBER_OF_PROCESSORS%`;
	}
	$threads =~ s/\s+$//;
	$threads = 0 unless $threads =~ /^\d+$/;
	unless ($threads){
		print "Cannot determine the number of CPUs. Do single threading.\n";
		$threads = 1;
	}
}

unless ($eqText){ # generate Entrez query text #
	if ($exUncultured){
		$eqText = "all [filter] NOT(environmental samples[filter] OR metagenomes[orgn])";
	}
	if (@searchTaxids){
		for (my $i=0; $i<scalar(@searchTaxids); $i++){
			next unless ($searchTaxids[$i] =~ /^\d+$/);
			$eqText .= " OR " if $eqText;
			$eqText .= "txid$searchTaxids[$i]\[orgn\]";
		}
	}
	if (@ignoreTaxids){
		$eqText .= " NOT (";
		for (my $i=0; $i<scalar(@ignoreTaxids); $i++){
			next unless $ignoreTaxids[$i] =~ /^\d+$/;
			$eqText .= " OR " if $i;
			$eqText .= "txid$ignoreTaxids[$i]\[orgn\]";
		}
		$eqText .= ")";
	}
}


## read query list ##

print "Reading input data...\n";
die "\nError: The input/ folder is not found.\n" unless (-d "$wkDir/input");
opendir (DIR, "$wkDir/input");
foreach my $file (grep {!/^\./} readdir DIR){
	next unless -s "$wkDir/input/$file";
	$ins{stem_name($file)} = { 'file'=>$file, 'prots'=>[] }
}
closedir DIR;
die "Error: No data are found in the input/ folder.\n" unless (%ins);
foreach my $set (sort keys %ins){
	$s = $ins{$set}{'file'};
	my $intype = ""; # type of input file format
	open IN, "<$wkDir/input/$s" or die "\nError: Failed to read input file $s.\n";
	while(<IN>){
		s/\s+$//; next unless $_; next if /^#/;
		unless ($intype){
			$intype = "list";
			$intype = "fasta" if (/^>/);
		}
		if ($intype eq 'fasta'){
			if (s/^>//){
				@a = seq_title ($_);
				push @{$ins{$set}{'prots'}}, { 'gi'=>$a[0], 'accn'=>$a[1], 'name'=>$a[2], 'product'=>$a[3], 'seq'=>'', 'hits'=>[], 'done'=>0 };
			}else{
				$ins{$set}{'prots'}[-1]{'seq'} .= $_; # append sequence
			}
		}elsif ($intype eq "list"){
			@a = seq_title ($_);
			push @{$ins{$set}{'prots'}}, { 'gi'=>$a[0], 'accn'=>$a[1], 'name'=>$a[2], 'product'=>$a[3], 'seq'=>'', 'hits'=>[], 'done'=>0 };
		}
	}
	close IN;
	my $i = scalar @{$ins{$set}{'prots'}};
	die "Error: sample $set does not contain any protein entries.\n" unless $i;
	print "  $set: $i proteins.\n";
}
$nProt += scalar @{$ins{$_}{'prots'}} for (keys %ins);
print "Done. $nProt proteins from ".(scalar keys %ins)." set(s) to query.\n";


## read search results from previous runs ##

if (-d "$wkDir/search"){
	print "Reading search results from previous run(s)...\n";
	foreach my $set (keys %ins){
		next unless -d "$wkDir/search/$set";
		my ($goodResults, $goodHits) = (0, 0);
		my ($badResults, $badHits) = (0, 0);
		foreach (my $i=0; $i<scalar(@{$ins{$set}{'prots'}}); $i++){
			my $file = "$wkDir/search/$set/".$ins{$set}{'prots'}[$i]{'name'}.".txt";
			if (-e $file){
				if (-s $file){
					my @out = ();
					open IN, "<$file" or next;
					while (<IN>){
						s/\s+$//;
						push (@out, $_);
					}
					close IN;
					my $badHitsHere = 0;
					my $k = 1; # status of report
					foreach (@out){
						next unless $_;
						if ($k == 1 and (/^#NEXUS/)){ $k = 2; next; }
						if ($k == 2 and /^BEGIN QUERY;/){ $k = 3; next; }
						if ($k == 3 and /^END;$/){ $k = 4; next; }
						if ($k == 4 and /^BEGIN ORGANISM;$/){ $k = 5; next; }
						if ($k == 5 and /^END;$/){ $k = 6; last; }
						if ($k == 5){
							next if /^\[/;
							next if /^;/;
							my @a = split (/\t/);
							if (@a and scalar(@a) < 6){
								$_ = "[deleted]";
								$badHitsHere ++;
								$badHits ++;
							}else{
								$goodHits ++;
							}
						}
					}
					if ($k < 6){ # bad report
						$badResults ++;
						unlink $file;
						print "  Warning: Incomplete search result: $set/".$ins{$set}{'prots'}[$i]{'name'}.". Deleted.\n";
					}else{
						$ins{$set}{'prots'}[$i]{'done'} = 1;
						$goodResults ++;
						if ($badHitsHere){
							open OUT, ">$file";
							foreach (@out){
								next if /^\[deleted\]/;
								print OUT "$_\n";
							}
							close OUT;
							print "  Warning: $badHitsHere invalid hits deleted from $set/".$ins{$set}{'prots'}[$i]{'name'}.".\n";
						}
					}
				}else{
					$badResults ++;
					unlink $file;
					print "  Empty search result: $set/".$ins{$set}{'prots'}[$i]{'name'}.". Deleted.\n";
				}
			}
		}
		print "  $set: $goodHits valid hits for $goodResults proteins.\n";
		$ins{$set}{'done'} = 1 if ($goodResults == scalar(@{$ins{$set}{'prots'}}));
		$nDone += $goodResults;
	}
	print "Done. $nDone results found, remaining ".($nProt - $nDone)." proteins to search.\n";
}

if ($nProt-$nDone <= 0){
	print "Batch homology search completed. searcher.pl exits.\n";
	print "You may proceed with HGT prediction by running analyzer.pl.\n\n";
	exit 0;
}


## detect pre-computed search results ##
# This function was inspired by Conor Meehan (cmeehan@itg.be) #

if ($preSearch){
	my $n = 0;
	opendir (DIR, $preSearch);
	foreach my $file (grep {!/^\./} readdir DIR){
		next unless -s "$preSearch/$file";
		$s = stem_name($file);
		if (exists $ins{$s}){
			$ins{$s}{'prefile'} = $file;
			$n ++;
		}
	}
	closedir DIR;
	print "Pre-computed search results are found for $n protein set(s).\n";
}

unless ($httpBlast){

	## read local taxonomy database ##
	
	unless ($protdb){ print "Warning: A protein database is required for local sequence homology search.\n"; }
	unless ($taxdump){ die "Error: A taxonomy database is required for local sequence homology search.\n"; }
	unless (-d $taxdump){ die "Error: Invalid taxonomy database directory: $taxdump.\n"; }
	unless (-s "$taxdump/nodes.dmp" and -s "$taxdump/names.dmp"){ die "Error: Taxonomy database is not found under $taxdump.\n" ; }
	print "Reading taxonomy database...";
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

	## read protein-to-TaxID dictionary ##
	
	if ($prot2taxid){
		unless (-s $prot2taxid){ die "Error: Invalid protein-to-TaxID dictionary: $prot2taxid.\n"; }
		print "Reading protein-to-TaxID dictionary...";
		open IN, "<$prot2taxid";
		while (<IN>){
			s/\s+$//; next unless $_;
			@a = split (/\s+/);
			next unless $#a;
			$prot2taxids{$a[0]} = $a[1];
		}
		close IN;
		print " done. ".scalar(keys %prot2taxids)." records read.\n";
	}elsif ($searchTool ne "BLAST"){
		print "Warning: A protein-to-TaxID dictionary is not provided.\n" ;
		if ($interactive){
			print "Press Enter to proceed, or Ctrl+C to exit:";
			$s = <STDIN>;
		}
	}
}


## read taxonomic information ##

if (-d "$wkDir/taxonomy"){
	print "Reading taxonomy records from previous run(s)...";
	if (-s "$wkDir/taxonomy/taxa.db"){
		open IN, "<$wkDir/taxonomy/taxa.db";
		while (<IN>){
			s/\s+$//; next if /^#/; next unless $_;
			@a = split (/\t/);
			next if exists $dbTaxa{$a[0]};
			my %taxon :shared = ('organism'=>$a[1], 'lineage'=>$a[2]);
			my $i = 3; $taxon{$_} = $a[$i++] for (@ranks);
			$dbTaxa{$a[0]} = \%taxon;
		}
		close IN;
	}
	if (-s "$wkDir/taxonomy/ranks.db"){
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

if (-s "$wkDir/taxonomy/self.info"){
	# The file self.info is like: sample_name | TaxID | organism_name
	my $nRead = 0;
	open IN, "<$wkDir/taxonomy/self.info";
	while (<IN>){
		s/\s+$//; next if /^#/; next unless $_;
		@a = split (/\t/);
		if (exists $ins{$a[0]}){
			$ins{$a[0]}{'taxid'} = $a[1];
			$ins{$a[0]}{'organism'} = $a[2];
		}
	}
	close IN;
}
if ($selfTax){ # TaxIDs of input protein sets
	%h = ();
	@a = split (/,/, $selfTax);
	foreach (@a){
		@b = split (/:/);
		next unless $#b;
		$h{$b[0]} = $b[1];
	}
	foreach my $set (keys %ins){
		next if exists $ins{$set}{'taxid'};
		next unless exists $h{$set};
		@a = get_taxonomy (($h{$set}));
		$ins{$set}{'taxid'} = $h{$set};
		$ins{$set}{'organism'} = $a[0];
		open OUT, ">>$wkDir/taxonomy/self.info";
		print OUT $set,"\t",$ins{$set}{'taxid'},"\t",$ins{$set}{'organism'},"\n";
		close OUT;
	}
}
if ($interactive){
	foreach my $set (keys %ins){
		next if exists $ins{$set}{'taxid'};
		print "Enter the TaxID of $set, or press Enter if you don't know:";
		$s = <STDIN>; chomp $s; next unless $s;
		@a = get_taxonomy (($s));
		$ins{$set}{'taxid'} = $s;
		$ins{$set}{'organism'} = $a[0];
		open OUT, ">>$wkDir/taxonomy/self.info";
		print OUT $set,"\t",$ins{$set}{'taxid'},"\t",$ins{$set}{'organism'},"\n";
		close OUT;
	}
}
print "Taxonomy of input protein sets:\n";
my @missingTax = ();
foreach my $set (sort keys %ins){
	if (exists $ins{$set}{'taxid'}){
		print "  $set: ", $ins{$set}{'organism'}, " (", $ins{$set}{'taxid'}, ")\n";
	}else{
		push (@missingTax, $set);
	}
}
if (@missingTax){
	print "Attempting to identify taxonomy of ".(scalar @missingTax)." protein set(s) :\n";
	print "  ", join (",", @missingTax), "\n";
	foreach my $set (@missingTax){
		my $taxid = "";
		if (%prot2taxids){ # look up the dictionary
			foreach my $prot (@{$ins{$set}{'prots'}}){
				$taxid = '';
				my ($gi, $accn, $name) = ($prot->{'gi'}, $prot->{'accn'}, $prot->{'name'});
				if (exists $prot2taxids{$name} and $prot2taxids{$name} !~ /,/){ $taxid = $prot2taxids{$name}; }
				elsif (exists $prot2taxids{$gi} and $prot2taxids{$name} !~ /,/){ $taxid = $prot2taxids{$gi}; }
				elsif (exists $prot2taxids{$accn} and $prot2taxids{$name} !~ /,/){ $taxid = $prot2taxids{$accn}; }
				next unless $taxid;
				next if $taxid =~ /,/;
				last;
			}
		}
		unless ($taxid){ # look up the BLAST database
			my ($gi, $accn, $name) = ($ins{$set}{'prots'}[0]{'gi'}, $ins{$set}{'prots'}[0]{'accn'}, $ins{$set}{'prots'}[0]{'name'});
			if ($searchTool eq "BLAST" and not $httpBlast){
				my $query = "";
				
				if ($accn){ $query = $accn; }
				elsif ($gi){ $query = $gi; }
				else { $query = $name; }
				my @out = `$blastdbcmd -db $protdb -entry $query -outfmt \"%a %g %T %t\"`;

				# The command is to query the BLAST database for one or more particular sequences.
				# The four codes represent accn, gi, taxid, title
				# The query may be gi or accn or "gi|###|ref|###" or "title".
				# If the database does not contain TaxID (or if not -taxid_map), %T will be 0.
				# when making database (using makeblastdb), one should do -parse_seqids to enable gi and accn search

				my $found = 1;
				foreach (@out){ if (/not found in BLAST database/){ $found = 0; last; } }
				if ($found){
					@a = split (/\s+/, $out[0]);
					$taxid = $a[2] if (scalar(@a) > 3 and $a[2] and $a[2] =~ /^\d+$/);
				}
				unless ($taxid){
					if ($gi){ # look up the NCBI server
						my $iRetry = 0;
						sleep 1;
						while (1){
							$s = get "$eSummaryServer?db=protein&id=$gi";
							last if (defined $s) and ($s =~ /<eSummaryResult>/s);
							die "\nFailed to retrieve taxonomic information from NCBI.\n" if ($iRetry >= $retries);
							$iRetry ++; sleep $delay; next;
						}
						$taxid = $1 if ($s =~ /<Item Name=\"TaxId\" Type=\"Integer\">(\d+?)<\/Item>/);
					}
				}
			}
		}
		if ($taxid){
			@a = get_taxonomy(($taxid));
			$ins{$set}{'taxid'} = $taxid;
			$ins{$set}{'organism'} = $a[0];
			open OUT, ">>$wkDir/taxonomy/self.info";
			print OUT $set, "\t", $taxid, "\t", $a[0], "\n";
			close OUT;
			print OUT "  $set: $taxid ($a[0])\n";
		}else{
			die "Error: Cannot identify the taxonomy of $set.\n";
		}
	}
}


## perform batch sequence similarity search ##

mkdir "$wkDir/search" unless -d "$wkDir/search";
if ($requests == 1 or not $httpBlast){
	foreach my $set (sort keys %ins){
		next if exists $ins{$set}{'done'};
		print "Batch homology search of $set (".(scalar @{$ins{$set}{'prots'}})." queries) started.\n";
		mkdir "$wkDir/search/$set" unless -d "$wkDir/search/$set";
		@a = localtime(time); $s = (++$a[4])."/$a[3]/".($a[5]+=1900)." $a[2]:$a[1]:$a[0]";
		open LOG, ">>$wkDir/search/$set.log";
		print LOG "Program started at $s.\n";
		print LOG "Number of queries: ". (scalar @{$ins{$set}{'prots'}}) .".\n";
		close LOG;
		
		unless ($httpBlast){ # local search
		
			%h = ();
			for (my $i=0; $i<scalar(@{$ins{$set}{'prots'}}); $i++){
				$h{$ins{$set}{'prots'}[$i]{'name'}} = $i unless $ins{$set}{'prots'}[$i]{'seq'};
			}
			if (%h){
				if ($searchTool ne "BLAST"){ die "Error: Cannot extract query sequences from a non-BLAST local database.\n"; }
				open OUT, ">$wkDir/seqids.txt"; 
				print OUT $_."\n" for (keys %h);
				close OUT;
				my @out = `$blastdbcmd -dbtype=prot -db $protdb -entry_batch $wkDir/seqids.txt -outfmt \"%a %s\"`;
				unlink "$wkDir/seqids.txt";
				if (join('', @out) =~ /not found in BLAST database/){ die "Error: Query sequences not found in $protdb.\n"; }
				foreach (@out){
					s/\s+$//; next unless $_;
					next if /^Error/;
					@a = split (/\s+/); next unless $#a;
					$a[0] =~ s/\.\d+$//;
					if (exists $h{$a[0]}){ $ins{$set}{'prots'}[$h{$a[0]}]{'seq'} = $a[1]; }
				}
			}
			%h = ();
			for (my $i=0; $i<scalar(@{$ins{$set}{'prots'}}); $i++){
				$h{$ins{$set}{'prots'}[$i]{'name'}} = $i unless $ins{$set}{'prots'}[$i]{'seq'};
			}
			if (%h){ die "Error: One or more query sequences cannot be extracted from $protdb.\n"; }

			# new feature: batch self search #
			
			if ($searchTool eq "RAPSearch2" or $searchTool eq "DIAMOND"){
				%h = ();
				open OUT, ">$wkDir/tmp.in";
				for (my $i=0; $i<scalar(@{$ins{$set}{'prots'}}); $i++){
					$h{$ins{$set}{'prots'}[$i]{'name'}} = $i;
					print OUT ">".$ins{$set}{'prots'}[$i]{'name'}."\n".$ins{$set}{'prots'}[$i]{'seq'}."\n";
				}
				close OUT;
				if ($searchTool eq "RAPSearch2"){
					`$prerapsearch -d $wkDir/tmp.in -n $wkDir/raptmp`;
					`$rapsearch -q $wkDir/tmp.in -d $wkDir/raptmp -t a -s f -o $wkDir/tmp -z $threads`;
					unlink "$wkDir/tmp.in";
					unlink "$wkDir/raptmp";
					unlink "$wkDir/raptmp.info";
					die "Error in running RAPSearch2. Please check." unless -s "$wkDir/tmp.m8";
					open IN, "<$wkDir/tmp.m8";
					while (<IN>){
						s/\s+$//; next unless $_; next if /^#/;
						@a = split (/\t/);
						next if ($#a < 11);
						next unless $a[0] eq $a[1];
						next unless exists $h{$a[0]};
						next if exists $ins{$set}{'prots'}[$h{$a[0]}]{'selfalign'};
						$ins{$set}{'prots'}[$h{$a[0]}]{'selfalign'} = $a[10]."/".$a[11]."/".$a[2]."/".length($ins{$set}{'prots'}[$h{$a[0]}]{'seq'});
					}
					close IN;
					unlink "$wkDir/tmp.aln";
					unlink "$wkDir/tmp.m8";
				}elsif ($searchTool eq "DIAMOND"){
					`$diamond makedb --in $wkDir/tmp.in -d $wkDir/tmp`;
					`$diamond blastp -p $threads -q $wkDir/tmp.in -d $wkDir/tmp -a $wkDir/tmp -t $wkDir`;
					my @out = `$diamond view -a $wkDir/tmp.daa`;
					unlink "$wkDir/tmp.in";
					unlink "$wkDir/tmp.dmnd";
					unlink "$wkDir/tmp.daa";
					die "Error in running DIAMOND. Please check." unless @out;
					foreach (@out){
						s/\s+$//;
						@a = split (/\t/);
						next unless $#a == 11;
						next unless $a[0] eq $a[1];
						next unless exists $h{$a[0]};
						next if exists $ins{$set}{'prots'}[$h{$a[0]}]{'selfalign'};
						$ins{$set}{'prots'}[$h{$a[0]}]{'selfalign'} = $a[10]."/".$a[11]."/".$a[2]."/".length($ins{$set}{'prots'}[$h{$a[0]}]{'seq'});
					}
				}
			}
		
			# end of new feature #
			
			my @ids = ();
			for (my $i=0; $i<scalar(@{$ins{$set}{'prots'}}); $i++){
				next if $ins{$set}{'prots'}[$i]{'done'};
				push (@ids, $i);
				next if exists $ins{$set}{'prefile'};
				next unless $queries;
				unless (scalar(@ids) % $queries){
					local_search ($set, \@ids);
					@ids = ();
				}
			}
			local_search ($set, \@ids) if (@ids);
			
		}else{ # single-process http BLAST #
		
			for (my $i=0; $i<scalar(@{$ins{$set}{'prots'}}); $i++){
				my $query = $ins{$set}{'prots'}[$i]{'name'};
				next if -s "$wkDir/search/$set/$query.txt";
				my $return = 0;
				print "  BLASTing $query... ";
				my $iRetry = 0;
				while (1){
					$return = http_blast ($set, $i);
					last if $return =~ /\/1$/;
					print "failed.\n";
					print "\n" and last if ($iRetry >= $retries);
					print "  Retrying...";
					$iRetry ++; sleep $delay;
				}
				open LOG, ">>$wkDir/search/$set.log";
				if ($return =~ /\/(\d+)\/1$/){
					print LOG "$query\t$1\n";
					print "done. $1 hits.\n";
				}
				else{ print LOG "$query\tfailed\n"; }
				close LOG;
				sleep $delay;
			}
		}
		@a = localtime(time); $s = (++$a[4])."/$a[3]/".($a[5]+=1900)." $a[2]:$a[1]:$a[0]";
		open LOG, ">>$wkDir/search/$set.log";
		print LOG "Program ended at $s.\n";
		close LOG;
		unlink "$wkDir/tmp.in";
		print "Batch homology search of $set (".(scalar @{$ins{$set}{'prots'}})." queries) completed.\n";
	}
	
}else{ # multi-process http BLAST #

	print "$requests http BLAST threads are running in parallel.\n";
	
	my %running = (); # accn -> 1
	my %retry = (); # accn -> times of retry
	my %rSets = (); # protein sets that have been started
	my %failed = ();
	
	while ($nDone + (scalar keys %failed) < $nProt){
		if (my $n = $requests-scalar(threads->list(threads::running))){
			for (1..$n){
				foreach my $set (keys %ins){
					next if $ins{$set}{'done'};
					unless (exists $rSets{$set}){
						print "Batch BLAST of $set (".(scalar @{$ins{$set}{'prots'}})." queries) started.\n";
						$rSets{$set} = 1;
						mkdir "$wkDir/search/$set" unless -d "$wkDir/search/$set";
						@a = localtime(time); $s = (++$a[4])."/$a[3]/".($a[5]+=1900)." $a[2]:$a[1]:$a[0]";
						open LOG, ">>$wkDir/search/$set.log"; print LOG "Program started at $s.\nNumber of queries: ". (scalar @{$ins{$set}{'prots'}}) .".\n"; close LOG;
					}
					my $started = 0;
					for (my $i=0; $i<scalar(@{$ins{$set}{'prots'}}); $i++){
						my $query = $ins{$set}{'prots'}[$i]{'name'};
						next if -s "$wkDir/search/$set/$query.txt";
						next if exists $running{$query};
						next if exists $failed{"$set|$query"};
						my $thr = threads->create(\&http_blast, $set, $i);
						print "  BLAST of $query started.\n";
						$running{$query} = 1;
						$started = 1;
						select (undef, undef, undef, 0.25); # this is to delay 0.25 seconds
						last;
					}
					last if $started;
					print "Batch BLAST of $set (".(scalar @{$ins{$set}{'prots'}})." queries) completed.\n";
					$ins{$set}{'done'} = 1;
					delete $rSets{$set};
					@a = localtime(time); $s = (++$a[4])."/$a[3]/".($a[5]+=1900)." $a[2]:$a[1]:$a[0]";
					open LOG, ">>$wkDir/search/$set.log"; print LOG "Program ended at $s.\n"; close LOG;
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
					open LOG, ">>$wkDir/search/$set.log";
					print LOG "$query\t$hits\n";
					close LOG;
				}else{
					$retry{$query} = 0 unless exists $retry{$query};
					$retry{$query} ++;
					if ($retry{$query} >= $retries){
						delete $retry{$query};
						$failed{"$set|$query"} = 1;
						print "  BLAST of $query failed.\n";
						open LOG, ">>$wkDir/search/$set.log"; print LOG "$query\tfailed\n"; close LOG;
					}else{
						print "  BLAST of $query failed. Scheduled to retry.\n";
					}
				}
			}
		}
		sleep $delay;
		if ($retryFailed and %failed and ($nDone + (scalar keys %failed) == $nProt)){
			$retryFailed = 0;
			%failed = ();
			$ins{$_}{'done'} = 0 for (keys %ins);
		}
	}
}
print "Batch homology search completed. searcher.pl exits.\n";
print "You may re-run searcher.pl to validate the results and finish incomplete searches.\n";
print "Or you may proceed with HGT prediction by running analyzer.pl.\n\n";
exit 0;


## sequence homology search using a local program ##
  # usage: local_search ($set, @ids)
  # the function will extract $ins{$set}{'prots'}[$id]'s for the search
  # if the pre-computed results are present, the function will read instead of search
  # return: 0 (success) or 1 (fail)

sub local_search {
	my ($set, $refID) = @_;
	my @ids = @$refID;
	my %name2id = (); # protein name (accn) to index
	my $outfile = $wkDir."/tmp.out"; # search result file name
	if (exists $ins{$set}{'prefile'}){ # read pre-computed results
		print "  Importing pre-computed search results of $set...";
		$outfile = $preSearch."/".$ins{$set}{'prefile'};
		$name2id{$ins{$set}{'prots'}[$_]{'name'}} = $_ for (@ids);
	}else{ # de novo search
		if (scalar(@ids) == scalar(@{$ins{$set}{'prots'}})){
			print "  $searchTool"."ing all proteins of $set...";
		}else{
			if (scalar(@ids) <= 32){
				@a = ();
				push (@a, $ins{$set}{'prots'}[$_]{'name'}) for (@ids);
				print "  $searchTool"."ing ".join(",", @a)."...";
			}else{
				@a = ();
				for (my $i=0; $i<3; $i++){
					push (@a, $ins{$set}{'prots'}[$ids[$i]]{'name'});
				}
				print "  $searchTool"."ing ".join(",", @a)."... (".scalar(@ids)." proteins)";
			}
		}
		unlink "$wkDir/tmp.in" if -e "$wkDir/tmp.in";
		open OUT, ">$wkDir/tmp.in";
		foreach my $id (@ids){
			print OUT ">".$ins{$set}{'prots'}[$id]{'name'}."\n".$ins{$set}{'prots'}[$id]{'seq'}."\n";
			$name2id{$ins{$set}{'prots'}[$id]{'name'}} = $id;
		}
		close OUT;
		if ($searchTool eq "BLAST"){
			$s = "$blastp -query $wkDir/tmp.in -db $protdb -out $outfile";
			$s .= " -num_threads $threads" if ($threads > 1);
			$s .= " -evalue $evalue" if ($evalue);
			$s .= " -max_target_seqs $nHits" if ($nHits);
			$s .= " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq\"";
			`$s`;
			# this is the standard tabular format, plus aligned part of subject sequence
		}elsif ($searchTool eq "RAPSearch2"){
			$s = "$rapsearch -q $wkDir/tmp.in -d $protdb -o $wkDir/tmp -t a -s f";
			$s .= " -z $threads" if ($threads > 1);
			$s .= " -e $evalue" if ($evalue);
			$s .= " -v $nHits" if ($nHits);
			`$s`;
			$outfile = $wkDir."/tmp.m8";
			unlink $wkDir."/tmp.aln";
		}elsif ($searchTool eq "DIAMOND"){
			$s = "$diamond blastp -p $threads -q $wkDir/tmp.in -d $protdb -a $wkDir/tmp -t $wkDir";
			$s .= " -e $evalue" if ($evalue);
			$s .= " -k $nHits" if ($nHits);
			`$s`;
			$s = "$diamond view -a $wkDir/tmp.daa -o $outfile";
			`$s`;
			unlink $wkDir."/tmp.daa";
		}else{
			unlink "$wkDir/tmp.in";
			die "Error: Search tool not specified and pre-computed results not found for $set.\n";
		}
		unlink "$wkDir/tmp.in";
	}
	unless (-e $outfile){
		print "Warning: Search result missing.\n";
		return 1;
	}
	
	# read search result file
	my (%accn2gi, %gi2accn) = ((), ());
	my %seqid2taxid = (); # need to look up the TaxIDs of these sequence IDs (GI, accn or name)
	open IN, "<$outfile";
	while (<IN>){
		s/\s+$//; next unless $_; next if /^#/;
		@a = split (/\t/);
		next unless exists $name2id{$a[0]};
		my $id = $name2id{$a[0]};
		$a[1] =~ s/\s+$//; # RAPSearch2 adds ^M to subject ID. probably a bug.
		next if ($evalue and $a[10] ne "*" and $evalue < $a[10]); # evalue cutoff
		next if ($identity and $a[2] ne "*" and $identity > $a[2]); # % identity cutoff
		next if ($coverage and $a[6] ne "*" and $a[7] ne "*" and $coverage > ($a[7]-$a[6]+1)/length($ins{$set}{'prots'}[$id]{'seq'})*100); # % coverage cutoff
		%h = ('sseqid'=>$a[1], 'pident'=>$a[2], 'evalue'=>$a[10], 'bitscore'=>$a[11], 'qstart'=>$a[6], 'qend'=>$a[7], 'sseq'=>'', 'name'=>$a[1]);
		$h{'sseq'} = $a[12] if ($#a >= 12 and $getAln);
		if ($h{'sseqid'} =~ /^gi\|(\d+)\|.+\|([A-Z0-9_]+)\.\d+\|(.*)$/){
			$h{'name'} = $2;
			$accn2gi{$2} = $1 unless exists $accn2gi{$2};
			$gi2accn{$1} = $2 unless exists $gi2accn{$1};
			$seqid2taxid{$1} = '0' unless exists $seqid2taxid{$1};
			$seqid2taxid{$2} = '0' unless exists $seqid2taxid{$2};
		}else{
			$seqid2taxid{$h{'name'}} = '0' unless exists $seqid2taxid{$h{'name'}};
		}
		push (@{$ins{$set}{'prots'}[$id]{'hits'}}, {%h});
	}
	close IN;
	unlink $outfile unless exists $ins{$set}{'prefile'};
	
	# get TaxIDs for hits
	if (%seqid2taxid){ 
		my @seqids4db = (); # sequence IDs not in dictionary
		
		# look up in the dictionary
		if (%prot2taxids){
			foreach my $seqid (keys %seqid2taxid){
				next if $seqid2taxid{$seqid};
				if (exists $prot2taxids{$seqid}){
					$seqid2taxid{$seqid} = $prot2taxids{$seqid};
					if (exists $accn2gi{$seqid} and not $seqid2taxid{$accn2gi{$seqid}}){ $seqid2taxid{$accn2gi{$seqid}} = $prot2taxids{$seqid}; }
					if (exists $gi2accn{$seqid} and not $seqid2taxid{$gi2accn{$seqid}}){ $seqid2taxid{$gi2accn{$seqid}} = $prot2taxids{$seqid}; }
				}
			}
			foreach my $seqid (keys %seqid2taxid){
				push (@seqids4db, $seqid) unless $seqid2taxid{$seqid};
			}
		}else{
			@seqids4db = keys %seqid2taxid;
		}
		
		# find the rest in database
		if (@seqids4db and $searchTool eq "BLAST" and not $httpBlast){
			open OUT, ">$wkDir/seqids.txt"; 
			print OUT $_."\n" for (keys %seqid2taxid);
			close OUT;
			my @out = `$blastdbcmd -dbtype=prot -db $protdb -entry_batch $wkDir/seqids.txt -outfmt \"%a %g %T\"`;
			unlink "$wkDir/seqids.txt";
			foreach (@out){
				s/\s+$//; next unless $_;
				@a = split (/\s+/); next if $#a < 2;
				next unless $a[2];
				$a[0] =~ s/\.\d+$//;
				if (exists $seqid2taxid{$a[0]} and not $seqid2taxid{$a[0]}){ $seqid2taxid{$a[0]} = $a[2]; }
				elsif ($a[1] ne "N/A" and exists $seqid2taxid{$a[1]} and not $seqid2taxid{$a[1]}){ $seqid2taxid{$a[1]} = $a[2]; }
			}
		}
	}
	
	# get complete taxonomy information
	my %taxid2orgn = ();
	foreach my $seqid (keys %seqid2taxid){
		next unless $seqid2taxid{$seqid};
		@a = split (/,/, $seqid2taxid{$seqid}); # one SeqID may correspond to multiple TaxIDs
		foreach (@a){
			next if exists $taxid2orgn{$_};
			$taxid2orgn{$_} = '';
		}
	}
	if (%taxid2orgn){
		@a = keys %taxid2orgn;
		@b = get_taxonomy (@a);
		for (my $i=0; $i<scalar(@a); $i++){
			$taxid2orgn{$a[$i]} = $b[$i];
		}
	}

	# ignore invalid organism names
	my @unwanted = ();
	foreach my $taxid (keys %taxid2orgn){
		my $organism = $taxid2orgn{$taxid};
		if ($organism eq 'na'){ push (@unwanted, $taxid); }
		elsif ($organism =~ /^Unresolved/){ push (@unwanted, $taxid); }
		if ($taxonUCASE and ($organism !~ /^\[?[A-Z]/)){ # a valid organism name may start with "["
			push (@unwanted, $taxid);
		}
		if (@ignoreTaxa){
			foreach (@ignoreTaxa){
				if ($organism =~ /$_/){ push (@unwanted, $taxid); last; }
			}
		}
	}
	delete $taxid2orgn{$_} for @unwanted;

	# process hits for each query
	foreach my $id (@ids){

		# extract valid hits from search result #
		my @hits;
		my $isQueryIn = 0; # whether query is among subjects
		my $name = $ins{$set}{'prots'}[$id]{'name'};

		if (@{$ins{$set}{'prots'}[$id]{'hits'}}){
			my %usedTaxids = ();
			my %usedSpecies = ();
			
			# sort by bit score
			my %id2score = ();
			for (my $i=0; $i<scalar(@{$ins{$set}{'prots'}[$id]{'hits'}}); $i++){
				$id2score{$i} = $ins{$set}{'prots'}[$id]{'hits'}[$i]{'bitscore'};
			}
			foreach my $i (sort { $id2score{$b} <=> $id2score{$a} } keys %id2score){
				%h = %{$ins{$set}{'prots'}[$id]{'hits'}[$i]};
				my $taxids = '';
				if (exists $seqid2taxid{$h{'name'}}){ $taxids = $seqid2taxid{$h{'name'}}; }
				elsif (exists $accn2gi{$h{'name'}} and exists $seqid2taxid{$accn2gi{$h{'name'}}}){ $taxids = $seqid2taxid{$accn2gi{$h{'name'}}}; }
				next unless $taxids;
				foreach my $taxid (split (/,/, $taxids)){
					next unless exists $taxid2orgn{$taxid};
					my $organism = $taxid2orgn{$taxid};
					if ($ins{$set}{'prots'}[$id]{'name'} eq $h{'name'}){ # self-align result is always retained
						$isQueryIn = 1;
					}else{ # remove redundancy from other hits
						next if $taxid2orgn{$taxid} eq 'na';
						next if exists $usedTaxids{$taxid};
						if ($ignoreSubspecies){
							next unless exists $dbTaxa{$taxid};
							next unless exists $dbTaxa{$taxid}{'species'};
							my $species = $dbTaxa{$taxid}{'species'};
							next unless $species;
							next if exists $usedSpecies{$species};
							$usedSpecies{$species} = $dbRanks{$species};
						}
						$usedTaxids{$taxid} = 1;
					}
					my %hit = ('accn'=>$h{'name'}, 'expect'=>$h{'evalue'}, 'score'=>$h{'bitscore'}, 'identity'=>$h{'pident'}, 'coverage'=>'*', 'taxid'=>$taxid, 'organism'=>$organism, 'sequence'=>$h{'sseq'});
					if ($h{'qstart'} ne "*" and $h{'qend'} ne "*"){ $hit{'coverage'} = sprintf("%.2f", ($h{'qend'}-$h{'qstart'}+1)/length($ins{$set}{'prots'}[$id]{'seq'})*100); }
					push (@hits, {%hit});
					last if (scalar @hits >= $nHits);
				}
			}
		}

		# perform self search, in case not targetted in previous steps #
		unless ($isQueryIn){
			if (exists $ins{$set}{'prots'}[$id]{'selfalign'}){
				@a = split(/\//, $ins{$set}{'prots'}[$id]{'selfalign'});
			}else{
				@a = self_align ($name, $ins{$set}{'prots'}[$id]{'seq'});
			}
			if (@a != (0,0,0,0)){
				my %hit = ('accn'=>$name, 'expect'=>$a[0], 'score'=>$a[1], 'identity'=>$a[2], 'coverage'=>'100.00', 'taxid'=>$ins{$set}{'taxid'}, 'organism'=>$ins{$set}{'organism'}, 'sequence'=>$ins{$set}{'prots'}[$id]{'seq'});
				unshift (@hits, {%hit});
				pop @hits if (scalar @hits > $nHits);
			}
		}

		# output result
		my ($ntax, $nchar) = (0, 0);
		open (OUT, ">>$wkDir/search/$set/$name.txt");
		print OUT "#NEXUS\nBEGIN QUERY;\n";
		if ($ins{$set}{'prots'}[$id]{'accn'}){
			print OUT "\tGI=".$ins{$set}{'prots'}[$id]{'gi'}.";\n\tAccession=".$ins{$set}{'prots'}[$id]{'accn'}.";\n";
		}else{
			print OUT "\tName=".$ins{$set}{'prots'}[$id]{'name'}.";\n";
		}
		print OUT "\tLength=".length($ins{$set}{'prots'}[$id]{'seq'}).";\n";
		if ($ins{$set}{'prots'}[$id]{'product'}){
			print OUT "\tProduct=".$ins{$set}{'prots'}[$id]{'product'}.";\n";
		}
		print OUT "\tOrganism=".$ins{$set}{'organism'}.";\n";
		print OUT "END;\n\n";
		print OUT "BEGIN ORGANISM;\n";
		print OUT "[Accession\tOrganism\tTaxID\tBit-score\tE-value\t\%Identity\t\%Coverage]\n";
		for (my $i=0; $i<scalar(@hits); $i++){
			next if exists $hits[$i]{'delete'};
			$ntax ++;
			$nchar = length($hits[$i]{'sequence'}) if ($hits[$i]{'sequence'} and $nchar < length($hits[$i]{'sequence'}));
			print OUT $hits[$i]{'accn'}."\t".$hits[$i]{'organism'}."\t".$hits[$i]{'taxid'}."\t".$hits[$i]{'score'}."\t".$hits[$i]{'expect'}."\t".$hits[$i]{'identity'}."\t".$hits[$i]{'coverage'}."\n";
			last if $i > $maxHits;
		}
		print OUT ";\nEND;\n\n";
		if ($getAln){
			print OUT "BEGIN DATA;\n";
			print OUT "\tDIMENSIONS NTAX=$ntax NCHAR=$nchar;\n\tFORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\n\tMATRIX\n";
			for (my $i=0; $i<scalar(@hits); $i++){
				next if exists $hits[$i]{'delete'};
				@a = split (/\//, $hits[$i]{'accn'});
				print OUT $a[0]."\t".$hits[$i]{'sequence'}."\n";
				last if $i > $maxHits;
			}
			print OUT ";\nEND;\n\n";
		}
		close OUT;
		open LOG, ">>$wkDir/search/$set.log";
		print LOG "$name\t$ntax\n";
		close LOG;
	}
	print " done.\n";
	return 0;
}


## BLAST via http connection to NCBI server ##
  # usage: http_blast ($set, $id)
  # return: (query, set, number_of_hits, whether_successful (0 - failed, 1- successful)) separated by "/"
  # refer to: ftp://ftp.ncbi.nlm.nih.gov/blast/documents/web_blast.pl
  # and: http://www.ncbi.nlm.nih.gov/BLAST/Doc/urlapi.html
  # note: the Perl thread join function has some problems, that why I used this inconvenient way.

sub http_blast{
	my ($set, $id) = ($_[0], $_[1]);

	# generate query information
	my %self = %{$ins{$set}{'prots'}[$id]}; # copy the whole record
	$self{'length'} = 0; # need to know length
	$self{'length'} = length($self{'seq'}) if $self{'seq'}; # ideally, sequence is available, otherwise see below
	my $query = $self{'name'};

	# send BLAST request #
	my $isError = 0;
	my $url = "$blastServer?CMD=Put&PROGRAM=blastp&DATABASE=$httpDb&FILTER=m S";
	$url .= "&EXPECT=$evalue" if $evalue;
	$url .= "&MAX_NUM_SEQ=$nHits" if ($nHits != 100);
	$url .= "&EQ_TEXT=$eqText" if $eqText;
	my $querykey = $self{'name'};
	$querykey = $self{'seq'} if $self{'seq'}; # user-defined query sequence
	$url .= "&QUERY=$querykey";
	$s = get $url;
	
	# $url = substr($args,0,2048) if (length($args) > 2048); # In some situations, the max size of a URL is 2048 bytes
	# @a = localtime(time); print "Post: $a[2]:$a[1]:$a[0], ";
	my $starttime = time;
	
	my $rid;
	# @a = localtime(time); print "Respond: $a[2]:$a[1]:$a[0], ";
	if ($s =~ /^    RID = (.*$)/m){ $rid = $1; print "RID=$rid "; }
	else{ return "$query/$set/0/0"; } ##########################
	if ($s =~ /\s\((\d+) letters\)/){ $self{'length'} = $1; }
	if ($s =~ /^    RTOE = (.*$)/m){
		my $rtoe = $1;
		# print "RTOE=$rtoe, ";
		if ($rtoe and $rtoe =~ /^\d+$/){ sleep $rtoe; }else{ sleep $delay; }
	}else{ return "$query/$set/0/0"; }
	while (1){
		sleep 1;
		$s = get "$blastServer?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
		if ($s =~ /\s+Status=([A-Z]+)/m){
			if ($1 eq "WAITING"){
				if (time - $starttime > $timeout){
					$isError = 1; last;
				}else{
					print ". "; sleep $delay; next;
				}
			}
			elsif ($1 eq "FAILED"){ $isError = 1; last; }
			elsif ($1 eq "UNKNOWN"){ $isError = 1; last; }
			elsif ($1 eq "READY"){ @a = localtime(time); last; }
			else{ $isError = 1; last; }
		}else{
			$isError = 1; last;
		}
		
		# if ($s =~ /\s+Status=WAITING/m){ next; }
		# if ($s =~ /\s+Status=FAILED/m){ $isError = 1; last; }
		# if ($s =~ /\s+Status=UNKNOWN/m){ $isError = 1; last; }
		# if ($s =~ /\s+Status=READY/m){
		# 	if ($s =~ /\s+ThereAreHits=yes/m){ last; }
		# 	else{ last; } # no hits;
		# }
	}
	if ($isError){ return "$query/$set/0/0"; }
	
	# retrieve tabular report #
	sleep 1;
	$url = "$blastServer?CMD=Get&ALIGNMENT_VIEW=Tabular&FORMAT_TYPE=Text&RID=$rid";
	$url .= "&MAX_NUM_SEQ=$nHits&DESCRIPTIONS=$nHits" if ($nHits != 100);
	$s = get $url;
	if (($s !~ /# blastp/) or ($s !~ /# Query:\s/)){ return "$query/$set/0/0"; }
	my %hits = ();
	my $hitid = 0;
	foreach (split(/\n/, $s)){
		if (/^# Query/ and not $self{'seq'}){
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
	my $i = 0; # count
	@a = (); # all results
	@b = (); # subset of GIs
	foreach (keys %hits){
		$i ++;
		push (@b, $_);
		if ($i == 190){ # in some situations, a URI should not exceed ~2000 characters, which is approximately 190 GIs.
			sleep 1;
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
		sleep 1;
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
	$hitid = 0;
	my ($taxid, $key) = (0, '');
	foreach ( sort {$hits{$a}{'id'} <=> $hits{$b}{'id'}} keys %hits){
		if ($hits{$_}{'id'} != $hitid){
			$hitid = $hits{$_}{'id'};
			$taxid = $hits{$_}{'taxid'};
			$key = $_;
		}else{
			if ($hits{$_}{'taxid'} == $taxid){
				$hits{$key}{'accn'} .= "/".$hits{$_}{'accn'};
				delete $hits{$_};
			}
		}
	}

	# merge duplicated hits (proteins with same bit score from one organism) #
	if ($mergeDuplicates){
		my $score = 0;
		($taxid, $key) = (0, '');
		foreach ( sort {($hits{$b}{'score'} <=> $hits{$a}{'score'}) or ($hits{$a}{'taxid'} <=> $hits{$b}{'taxid'})} keys %hits){
			if (($hits{$_}{'score'} != $score) or ($hits{$_}{'taxid'} != $taxid)){
				$score = $hits{$_}{'score'};
				$taxid = $hits{$_}{'taxid'};
				$key = $_;
			}else{
				$hits{$key}{'accn'} .= "/".$hits{$_}{'accn'};
				delete $hits{$_};
			}
		}
	}

	# reorder accession numbers #
	foreach (keys %hits){
		$hits{$_}{'accn'} = order_accns ($hits{$_}{'accn'}, $query);
	}

	# retrieve taxonomy information from NCBI server #
	$i = 0;
	@a = ();
	foreach (keys %hits){
		next if exists $dbTaxa{$hits{$_}{'taxid'}};
		push @a, $hits{$_}{'taxid'};
		$i ++;
		if ($i == 150){
			sleep 1;
			get_taxonomy @a;
			$i = 0; @a = (); sleep $delay;
		}
	}
	sleep 1;
	get_taxonomy @a if @a;

	# check whether BLAST result contains query itself #
	my $isQueryIn = 0;
	foreach (keys %hits){
		@a = split (/\//, $hits{$_}{'accn'});
		foreach $s (@a){ if ($query eq $s){ $isQueryIn = $hits{$_}{'id'}; last; }}
		if ($isQueryIn){ $self{'length'} = $hits{$_}{'length'} unless $self{'length'}; last; }
	}

	# find out length of query sequence, in case not in previous steps #
	unless ($self{'length'} or $self{'seq'}){
		sleep 1;
		while (1){
			$s = get "$eSearchServer?db=protein&term=$query";
			last if (defined $s);
			sleep $delay;
		}
		$s =~ /<Id>(\d+)<\/Id>/;
		sleep 1;
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

	# perform self search, in case not targetted in previous steps #
	unless ($isQueryIn){
		sleep 1;
		@a = self_align ($query, '');
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
		my $i = 0;
		foreach (sort {$hits{$b}{'score'} <=> $hits{$a}{'score'}} keys %hits){
			$i ++;
			delete $hits{$_} if ($i > $maxHits);
		}
	}

	# create output file #
	open (OUT, ">$wkDir/search/$set/$query.txt");

	print OUT "#NEXUS\nBEGIN QUERY;\n";
	if ($self{'accn'}){
		print OUT "\tGI=".$self{'gi'}.";\n\tAccession=".$self{'accn'}.";\n";
	}else{
		print OUT "\tName=".$self{'name'}.";\n";
	}
	print OUT "\tLength=".$self{'length'}.";\n";
	if ($self{'product'}){
		print OUT "\tProduct=".$self{'product'}.";\n";
	}
	print OUT "\tOrganism=".$ins{$set}{'organism'}.";\n";
	print OUT "END;\n\n";

	# retrieve taxonomy report (using TaxBLAST)
	if ($taxBlast){
		sleep 1;
		$s = get "$blastServer?CMD=Get&FORMAT_TYPE=HTML&FORMAT_OBJECT=TaxBlast&RID=$rid&ALIGNMENTS=100000";
		if (($s !~ /Tax BLAST Report/) or ($s !~ /Lineage Report/)){
			print OUT "BEGIN ERROR;\nRetrieval of taxonomy report failed.\n;\nEND;\n\n";
		}else{
			my $reading = 0;
			foreach (split(/\n/, $s)){
				if ($_ eq "<B><A NAME=lineage>Lineage Report</A></B><BR>"){
					print OUT "BEGIN LINEAGE;\n";
					$reading = 1; next;
				}
				if (($_ eq "</FONT></PRE><HR>") and $reading){
					print OUT ";\nEND;\n\n";
					$reading = 0; next;
				}
				if ($reading){
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
		my $i = 0; # count
		my $allinfo = "";
		@b = (); # subset of GIs
		foreach (keys %hits){
			$i ++;
			push (@b, $_);
			if ($i == 190){ # a URI should not exceed ~2000 characters, which is approximately 190 GIs.
				sleep 1;
				while (1){
					$s = get $eFetchServer."?db=protein&rettype=FASTA&id=".join (",", @b);
					last if (defined $s);
					sleep $delay;
				}
				$allinfo .= $s;
				$i = 0; @b = ();
			}
		}
		if (@b){
			sleep 1;
			while (1){
				$s = get $eFetchServer."?db=protein&rettype=FASTA&id=".join (",", @b);
				last if (defined $s);
				sleep $delay;
			}
			$allinfo .= $s;
		}
		$i = 0; # current GI
		foreach (split (/\n/, $allinfo)){
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
		sleep 1;
		$s = get "$blastServer?CMD=Get&ALIGNMENT_VIEW=FlatQueryAnchoredNoIdentities&FORMAT_TYPE=Text&RID=$rid&ALIGNMENTS=100000";
		if (($s !~ /blastp/i) or ($s !~ /\nQuery=\s/)){
			print OUT "BEGIN ERROR;\nRetrieval of multiple sequence alignment failed.\n;\nEND;\n\n";
		}else{
			@c = split(/\n/, $s);
			my $reading = 0; # reading status
			my $iBlock = -1; # block ID
			my %seqs; # sequence alignment
			foreach (@c){
				if ($_ eq "ALIGNMENTS"){ $reading = 1; next; }
				if ($reading and /^\s/){ $reading = 0; last; }
				next unless $reading;
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
			foreach my $hit (keys %hits){
				@a = split /\//, $hits{$hit}{'accn'};
				foreach (@a){
					if (exists $seqs{$_}){
						$hits{$hit}{'sequence'} = $seqs{$_};
						last;
					}
				}
			}
		}
	}

	# output sequences #
	if ($seqBlast or $alnBlast){
		my $ntax = 0; # number of sequences
		my $nchar = 0; # maximum length of sequence
		foreach (keys %hits){
			next unless exists $hits{$_}{'sequence'};
			$ntax ++;
			$nchar = length ($hits{$_}{'sequence'}) if (length ($hits{$_}{'sequence'}) > $nchar);
		}
		print OUT "BEGIN DATA;\n";
		print OUT "\tDIMENSIONS NTAX=$ntax NCHAR=$nchar;\n\tFORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\n\tMATRIX\n";
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


## align the query sequence to itself ##
  # parameter: name, seq
  # return: (e-value, bit-score, identity, length)
  # return (0, 0, 0, 0) if failed

sub self_align {
	my ($name, $seq, $length) = ($_[0], $_[1], 0);
	$length = length($seq) if $seq;
	unless ($httpBlast){ # local mode
		my $fail = 0;
		my @result = ();
		open OUT, ">$wkDir/tmp.in";
		print OUT ">$name\n".$seq."\n";
		close OUT;
		if ($searchTool eq "BLAST"){
			my @out = `$blastp -query $wkDir/tmp.in -subject $wkDir/tmp.in -outfmt \"6 evalue bitscore pident\"`;
			unlink "$wkDir/tmp.in";
			@a = split (/\t/, $out[0]);
			if ($#a == 2){
				$a[2] =~ s/\s+$//;
				@result = ($a[0], $a[1], $a[2], $length);
			}else{ $fail = 1; }
		}elsif ($searchTool eq "RAPSearch2"){
			`$prerapsearch -d $wkDir/tmp.in -n $wkDir/raptmp`;
			`$rapsearch -q $wkDir/tmp.in -d $wkDir/raptmp -t a -s f -o $wkDir/tmp`;
			unlink "$wkDir/tmp.in";
			unlink "$wkDir/raptmp";
			unlink "$wkDir/raptmp.info";
			if (-s "$wkDir/tmp.m8"){
				open IN, "<$wkDir/tmp.m8";
				while (<IN>){
					s/\s+$//; next unless $_; next if /^#/;
					@a = split (/\t/);
					if ($#a < 11){ $fail = 1; last; }
					@result = ($a[10], $a[11], $a[2], $length);
					last;
				}
				close IN;
				unlink "$wkDir/tmp.aln";
				unlink "$wkDir/tmp.m8";
			}else{ $fail = 1; }

		}elsif ($searchTool eq "DIAMOND"){
			`$diamond makedb --in $wkDir/tmp.in -d $wkDir/tmp`;
			`$diamond blastp -p 1 -q $wkDir/tmp.in -d $wkDir/tmp -a $wkDir/tmp -t $wkDir`;
			my @out = `$diamond view -a $wkDir/tmp.daa`;
			unlink "$wkDir/tmp.in";
			unlink "$wkDir/tmp.dmnd";
			unlink "$wkDir/tmp.daa";
			if (@out){
				@a = split (/\t/, $out[0]);
				if ($#a == 11){
					$a[11] =~ s/\s+$//;
					@result = ($a[10], $a[11], $a[2], $length);
				}else{ $fail = 1; }
			}else{ $fail = 1; }
		}else{ $fail = 1; }
		if ($fail){
			unlink "$wkDir/tmp.in" if -e "$wkDir/tmp.in";
			return (0,0,0,0);
		}else{ return @result; }

	}else{ # http mode
		my $isError = 0;
		if ($seq){ $s = get "$blastAlignServer?CMD=Put&PROGRAM=blastp&DATABASE=$protdb&QUERY=".$seq."&SUBJECTS=".$seq; }
		else{ $s = get "$blastAlignServer?CMD=Put&PROGRAM=blastp&DATABASE=$protdb&QUERY=$name&SUBJECTS=$name"; };
		return (0,0,0,0) unless (defined $s and $s =~ /^    RID = (.*$)/m);
		my $rid = $1;
		if ($s =~ /\s\((\d+) letters\)/ and not $length){ $length = $1; }
		sleep 1;
		while (1){
			$s = get "$blastServer?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
			return (0,0,0,0) unless defined $s and $s =~ /\s+Status=(.+)/m;
			if ($1 eq "WAITING"){ sleep $delay; next; }
			if ($1 eq "FAILED" or $1 eq "UNKNOWN"){ $isError = 1; last; }
			if ($1 eq "READY" and $s =~ /\s+ThereAreHits=yes/m){ last; }
			$isError = 1; last;
		}
		return (0,0,0,0) if $isError;
		sleep 1;
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


## retrieve complete taxonomy information ##
  # parameter: array of TaxIDs
  # return: array of organism names in same order
  # other: write global variable %dbTaxa and %dbRanks, write file taxa.db and ranks.db

sub get_taxonomy{
	my @organisms = ();
	my %taxa2w = (); # taxa to write
	my %ranks2w = (); # ranks to write
	my %hRanks = map { $_ => 1 } @ranks;

	unless ($httpBlast){ # local taxonomy database
		foreach my $taxid (@_){
			if (exists $dbTaxa{$taxid}){
				push (@organisms, $dbTaxa{$taxid}{'organism'});
				next;
			}
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
				push (@organisms, "na");
				unless (exists $badTaxids{$taxid}){
					print " Warning: Invalid TaxID: $taxid.\n";
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
		my %taxid2organism = ();
		foreach (@_){
			if (exists $dbTaxa{$_}){ $taxid2organism{$_} = $dbTaxa{$_}{'organism'}; }
			else{ push (@taxids2get, $_) }
		}
		sleep 1;
		while (1){
			$s = get "$eFetchServer?db=taxonomy&id=".join (",", @taxids2get);
			last if (defined $s);
			die "\nFailed to retrieve taxonomic information from NCBI.\n" if ($iRetry >= $retries);
			$iRetry ++; sleep $delay; next;
		}
		if ($s =~ /<ERROR>ID list is empty/){
			print "Warning: Invalid TaxIDs:", join (",", @taxids2get), ".\n";
		}else{
			$s =~ s/<TaxaSet>//;
			while ($s =~ s/\n<Taxon>\s+<TaxId>(\d+)<\/TaxId>\s+<ScientificName>(.+?)<\/ScientificName>(.+?)\n<\/Taxon>//s){
				my $id = $1; my $t = $3;
				$taxid2organism{$id} = $2;
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
		foreach (@_){
			if (exists $taxid2organism{$_}){
				push (@organisms, $taxid2organism{$_});
			}else{
				push (@organisms, "na");
				print "Warning: Invalid TaxID: $_.\n";
				unless (exists $badTaxids{$_}){
					$badTaxids{$_} = 1;
					open OUT, ">>$wkDir/taxonomy/invalid.taxids";
					print OUT $_."\n";
					close OUT;
				}
			}
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

## parse sequence title ##
  # GenBank-style sequence titles have gi, accn, product and name (accn)
  # plain titles have name and product

sub seq_title {
	my $title = $_[0];
	my ($gi, $accn, $name, $product) = ('', '', '', '');
	if ($title =~ /^gi\|(\d+)\|.+\|([A-Z0-9_]+)\.\d+\|(.*)$/){ # GenBank-style title
		($gi, $accn, $name) = ($1, $2, $2);
		if ($3){
			$product = $3;
			$product = $1 if ($product =~ /^(.+)\[/);
			$product =~ s/^\s+|\s+$//g;
		}
	}else{ # user-defined name
		$title =~ s/^\s+|\s+$//g;
		if ($title =~ /^(\S+)\s+(.+)$/){
			($name, $product) = ($1, $2);
		}else{
			$name = $title;
		}
	}
	return ($gi, $accn, $name, $product);
}

# reorder accession numbers #
sub order_accns {
	# order: NP_ > XP_ = YP_ = ZP_ = AP_ > everything else
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

# get stem file name #
sub stem_name {
	my $stem = $_[0];
	$stem =~ s/\.[^\.]+$//;
	return $stem;
}

