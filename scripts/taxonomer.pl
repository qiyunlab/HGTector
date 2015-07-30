#!/usr/bin/perl

use strict;
use warnings;
$| = 1;

print "

This script retrieves hierarchical taxonomy information for BLAST hits
from a remote or local taxonomy database.

Usage:
  perl taxonomer.pl <working directory>

Output:
  taxonomy/self.info
  taxonomy/taxa.db
  taxonomy/ranks.db

" and exit unless @ARGV;

print "
-> Taxonomer: Retrieve taxonomic information for BLAST hits. <-
";

use LWP::Simple;

sub get_taxonomy;


## global variables ##

my $i; my $s; my $t; my @a; my %h;

my $wkDir = $ARGV[0];					# working directory
my %taxids = ();
my %dbTaxa = ();						# taxa.db
my %dbRanks = ();						# ranks.db
my %taxdumps = ();						# local taxonomy database: taxid => (name, parent, rank)


## program parameters ##

my $interactive = 1;					# interactive or automatic mode
my @ranks = ('species', 'genus', 'family', 'order', 'class', 'phylum');
my $blastMode = 0;						# BLAST mode (0: via http connection, 1: local BLAST)
my $nRetry = 5;							# maximum number of retries
my $delay = 5;							# time delay (seconds) between two http requests
my $eSearchServer = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
my $eFetchServer = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
my $eSummaryServer = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi";
my $taxdump = "";						# directory of NCBI taxonomy database (nodes.dmp and names.dmp)


## read configuration ##

if (-e "$wkDir/config.txt"){
	open IN, "<$wkDir/config.txt";
	while (<IN>){
		s/#.*$//; s/\s+$//g; s/^\s+//g; next unless $_;
		$interactive = $1 if /^interactive=([01])$/;
		@ranks = split (/,/, $1) if /^ranks=(.+)$/;
		$blastMode = $1 if /^blastMode=([01])$/;
		$nRetry = $1 if /^nRetry=(\d+)$/;
		$delay = $1 if /^delay=(\d+)$/;
		$eSearchServer = $1 if /^eSearchServer=(.+)$/;
		$eFetchServer = $1 if /^eFetchServer=(.+)$/;
		$eSummaryServer = $1 if /^eSummaryServer=(.+)$/;
		$taxdump = $1 if /^taxdump=(.+)$/;
	}
}


## overwrite old result ##

mkdir "$wkDir/taxonomy" unless -d "$wkDir/taxonomy";

if (-e "$wkDir/taxonomy/taxa.db" and -e "$wkDir/taxonomy/ranks.db"){
	print "Warning: taxa.db and ranks.db already exist. ";
	if ($interactive){
		print "Press Enter to overwrite, or Ctrl+C to exit:";
		$s = <STDIN>;
	}else{
		print "Overwritten.\n";
	}
	unlink "$wkDir/taxonomy/taxa.db";
	unlink "$wkDir/taxonomy/ranks.db";
}


## read protein sets ##

opendir (DIR, "$wkDir/blast");
my @sets = grep {-d "$wkDir/blast/$_" and not /^\.{1,2}$/} readdir(DIR);
close DIR;
die "No BLAST results detected.\n" unless @sets;
print "Found ".(scalar @sets)." sets of BLAST reports. Reading TaxIDs...\n";


## read self information ##

my %ins = map { $_ => 0 } @sets;
if (-e "$wkDir/taxonomy/self.info"){
	open IN, "<$wkDir/taxonomy/self.info";
	while (<IN>){
		s/\s+$//; next if /^#/; next unless $_;
		@a = split (/\t/); next unless @a;
		foreach (keys %ins){
			if ($a[0] eq $_){
				$ins{$_} = 1;
				$taxids{$a[2]} = 1 unless $a[2];
				last;
			}
		}
	}
	close IN;
}

my @set2get = (); my @taxon4set = ();
foreach (keys %ins){
	push (@set2get, $_) unless $ins{$_};
}

if ($interactive and @set2get){
	print "Enter the NCBI TaxIDs in order for the following protein sets, separated by space. Type 0 if not applicable:\n ";
	print " $_" for (@set2get);
	print "\n";
	$s = <STDIN>; chomp $s;
	if ($s){
		@taxon4set = split (/\s+/, $s);
		$taxids{$_} = 1 for (@taxon4set);
	}
}

foreach my $set (@sets){
	print " $set";
	opendir DIR, "$wkDir/blast/$set";
	my @blas = grep (s/\.bla$//, readdir(DIR));
	close DIR;
	foreach my $bla (@blas){
		my $reading = 0;
		open IN, "<$wkDir/blast/$set/$bla.bla";
		while (<IN>){
			s/\s+$//; next unless $_;
			next if /^#/; next if /^;/; next if /^\[/;
			if (/^BEGIN ORGANISM/){ $reading = 1; next; }
			if (/^END;/ and $reading){ last; }
			if ($reading){
				@a = split (/\t/);
				next if $#a < 5;
				next unless $a[2] =~ /^\d+$/;
				next unless $a[2];
				next if exists $taxids{$a[2]};
				$taxids{$a[2]} = 1;
			}
		}
	}
}
my $n = scalar (keys %taxids);
print "\nDone. $n TaxIDs read.\n";


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
	print " done.\n";
}


## retrieve taxonomy information ##

print "Retrieving taxonomy information...\n";
print "0-------------25-------------50------------75------------100%\n";

my ($count, $done, $progress) = (0, 0, 0); 
@a = ();
foreach my $taxid (sort {$a<=>$b} keys %taxids){
	push (@a, $taxid);
	$count ++; $done ++;
	if ($count == 150){
		get_taxonomy @a;
		if ($done/$n >= $progress/60){ 
			print "."; $progress ++;
		}
		@a = (); $count = 0; sleep $delay;
	}
}
get_taxonomy @a if @a;
print "done.\n";

if (@set2get and @taxon4set){
	open OUT, ">>$wkDir/taxonomy/self.info";
	for ($i=0; $i<=$#set2get; $i++){
		if ($#taxon4set >= $i and $taxon4set[$i] and exists $dbTaxa{$taxon4set[$i]}){
			print OUT $set2get[$i]."\tna\t".$taxon4set[$i]."\t".$dbTaxa{$taxon4set[$i]}{'organism'}."\n";
		}
	}
	close OUT;
}

print "Taxonomy information has been retrieved and saved.\n";

exit 0;


## retrieve taxonomy information from remote or local database ##
  # this function differs from its counterpart in blaster.pl in that it does not support multithreading.

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
				print "TaxID $taxid is not found in the local taxonomy database.\n";
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

