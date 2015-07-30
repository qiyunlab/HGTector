#!/usr/bin/perl

use strict;
use warnings;
$| = 1;

print "

This script identifies invalid or redundant hits in BLAST reports based on taxonomy. It
does not delete them, instead, it labels them with an 'x' so that the subsequent steps
will ignore them.

Usage:
  perl filter.pl <working directory>

Output:
  Will be appended to each BLAST report file (*.bla).

" and exit unless @ARGV;

print "
-> Filterer: Label invalid or redundant BLAST hits. <-
";

## global variables ##

my $i; my $j; my $n; my $s; my $t; my @a; my @b; my @c; my %h;

## program parameters ##

my $wkDir = $ARGV[0];							# working directory

my @ignoreTaxa = ();							# ignore taxon names containing the following words
my $taxonUCASE = 0;								# ignore taxon names that do not start with a capital letter
my $ignoreParalogs = 1;							# ignore potential paralogs (hits with same taxon names but different bit scores)
my $ignoreSeqRepeats = 1;						# ignore repeated sequences (hits targetting different regions of the same protein)
my $ignoreSubspecies = 0;						# ignore more than one subspecies from the same species
my @ranks = ('species', 'genus', 'family', 'order', 'class', 'phylum');

## read configurations ##

if (-e "$wkDir/config.txt"){
	open IN, "<$wkDir/config.txt";
	my $readMonitor = 0; my $readMiddle = 0;
	while (<IN>){
		s/#.*$//; s/\s+$//g; s/^\s+//g; next unless $_;
		@ignoreTaxa = split(/,/, $1) if /^ignoreTaxa=(.+)$/;
		$taxonUCASE = $1 if /^taxonUCASE=([01])$/;
		$ignoreParalogs = $1 if /^ignoreParalogs=([01])$/;
		$ignoreSeqRepeats = $1 if /^ignoreSeqRepeats=([01])$/;
		$ignoreSubspecies = $1 if /^ignoreSubspecies=([01])$/;
		@ranks = split (/,/, $1) if /^ranks=(.+)$/;
	}
	close IN;
}

## read taxonomy information ##

my %dbTaxa = ();
my %dbRanks = ();
if ($ignoreSubspecies and -d "$wkDir/taxonomy"){
	print "Reading taxonomy database...";
	if (-e "$wkDir/taxonomy/taxa.db"){
		open IN, "<$wkDir/taxonomy/taxa.db";
		while (<IN>){
			s/\s+$//; next if /^#/; next unless $_;
			@a = split (/\t/);
			next if exists $dbTaxa{$a[0]};
			my %taxon = ('organism',$a[1],'lineage',$a[2]);
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
}

## read protein sets ##

my @sets = ();
print "\nReading protein sets...";
opendir (DIR, "$wkDir/blast") or die "Blast folder not exist.\n";
@a = readdir(DIR);
close DIR;

foreach (@a){
	next if (/^\./);
	push @sets, $_ if -d "$wkDir/blast/$_";
}
print "done. ";
die "No protein sets detected.\n" unless @sets;
print @sets." protein sets detected.\n";

print "\n0-------------25-------------50------------75------------100%\n";

foreach my $set (@sets){
	my $dir = "$wkDir/blast/$set";
	opendir (DIR, "$dir") or next;
	my @files = grep (/\.bla$/, readdir(DIR));
	close DIR;
	next unless @files;
	
	my $iProtein = 0;
	my $iProgress = 0;
	my $nProtein = $#files+1;
	print "$set has $nProtein proteins. Filtering...\n";
	
	foreach my $file (@files){
		$iProtein ++;
		my @hits = ();
		my %self = ();
		my @infile = ();
		open IN, "<$dir/$file" or next;
		while (<IN>){
			s/\s+$//;
			push (@infile, $_);
		}
		close IN;

		my $reading = 0;
		my ($hasCoverage, $hasDistance) = (0, 0);
		foreach (@infile){
			s/\s+$//;
			if (/^BEGIN QUERY/){ $reading = "query"; next; }
			if (/^BEGIN ORGANISM/){ $reading = "organism"; next; }
			if (/^END;/){ $reading = 0; next; }
			if ($reading eq "query"){
				$self{'accn'} = $1 if /^\tName=(.+);$/;
				if (/^\tAccession=(.+);$/){ $self{'accn'} = $1; $self{'accn'} =~ s/\.[\d]+$//; }
			}
			if ($reading eq "organism"){ # read organisms
				last unless exists $self{'accn'};
				next if /^;/;
				if (/^\[/){
					$hasCoverage = 1 if /Coverage/;
					$hasDistance = 1 if /Distance/;
					next;
				}
				@a = split (/\t/);
				my %hit = ();
				$hit{'accns'} = $a[0];
				$hit{'organism'} = $a[1];
				$hit{'taxid'} = $a[2];
				$hit{'score'} = $a[3];
				$hit{'evalue'} = $a[4];
				$hit{'identity'} = $a[5];
				$hit{'coverage'} = $a[6] if $hasCoverage;
				if ($#a >= 6+$hasCoverage and $a[6+$hasCoverage] ne 'x'){ $hit{'distance'} = $a[6+$hasCoverage]; }
				@a = split(/\//, $a[0]);
				$hit{'accn'} = $a[0];
				push @hits, {%hit};
				unless (exists $self{'taxid'}){ # identify self
					foreach (@a){
						if ($self{'accn'} eq $_){
							$self{'id'} = $#hits;
							$self{'organism'} = $hit{'organism'};
							$self{'taxid'} = $hit{'taxid'};
							$self{'score'} = $hit{'score'};
							$self{'evalue'} = $hit{'evalue'};
							$self{'identity'} = $hit{'identity'};
							$self{'coverage'} = $hit{'coverage'} if exists $hit{'coverage'};
							$self{'distance'} = $hit{'distance'} if exists $hit{'distance'};
							last;
						}
					}
				}
			}
		}

		## ignore lower lettered taxon names ##
		
		if ($taxonUCASE){
			for ($i=0; $i<=$#hits; $i++){
				if ($hits[$i]{'organism'} !~ /^[A-Z]/){
					$hits[$i]{'ignore'} = 1;
				}
			}
		}

		## ignore user-defined taxon names ##
		
		if (@ignoreTaxa){
			for ($i=0; $i<=$#hits; $i++){
				$j = 0;
				foreach (@ignoreTaxa){
					if ($hits[$i]{'organism'} =~ /$_/){ $j = 1; last; }
				}
				$hits[$i]{'ignore'} = 1 if $j;
			}
		}

		## ignore paralogs ##
		
		if ($ignoreParalogs){
			my %usedTaxids = ();
			for ($i=0; $i<=$#hits; $i++){
				next if exists $hits[$i]{'ignore'};
				if (exists $usedTaxids{$hits[$i]{'taxid'}}){
					$hits[$i]{'ignore'} = 1;
				}else{
					$usedTaxids{$hits[$i]{'taxid'}} = 1;
				}
			}
		}

		## ignore repeated sequences ##
		
		if ($ignoreSeqRepeats){
			for ($i=1; $i<=$#hits; $i++){
				next if exists $hits[$i]{'ignore'};
				if ($hits[$i]{'accn'} eq $hits[$i-1]{'accn'}){
					$hits[$i]{'ignore'} = 1;
				}
			}
		}

		## ignore subspecies ##
		
		if ($ignoreSubspecies){
			my %speciesDB = ();
			for ($i=0; $i<=$#hits; $i++){
				unless (exists $dbTaxa{$hits[$i]{'taxid'}} and exists $dbTaxa{$hits[$i]{'taxid'}}{'species'} and $dbTaxa{$hits[$i]{'taxid'}}{'species'}){
					$hits[$i]{'ignore'} = 1; next;
				}
				$s = $dbTaxa{$hits[$i]{'taxid'}}{'species'};
				if (exists $speciesDB{$s}){
					next if ($hits[$i]{'taxid'} == $self{'taxid'});
					$t = $self{'accn'}; next if (($hits[$i]{'accn'} eq $t) or ($hits[$i]{'accn'} =~ /^$t\//) or ($hits[$i]{'accn'} =~ /\/$t$/) or ($hits[$i]{'accn'} =~ /\/$t\//));
					$hits[$i]{'ignore'} = 1;
				}else{ $speciesDB{$s} = $dbRanks{$s}; }
			}
		}

		## write output ##

		open OUT, ">$dir/$file" or next;
		foreach (@infile){
			print OUT $_."\n";
			last if /^BEGIN ORGANISM/;
		}
		print OUT "[Accession\tOrganism\tTaxID\tBit-score\tE-value\t\%Identity\tIgnore]\n";
		for ($i=0; $i<=$#hits; $i++){
			print OUT $hits[$i]{'accns'}."\t".$hits[$i]{'organism'}."\t".$hits[$i]{'taxid'}."\t".$hits[$i]{'score'}."\t".$hits[$i]{'evalue'}."\t".$hits[$i]{'identity'}."\t";
			print OUT $hits[$i]{'coverage'}."\t" if exists $hits[$i]{'coverage'};
			print OUT $hits[$i]{'distance'}."\t" if exists $hits[$i]{'distance'};
			print OUT "x" if exists $hits[$i]{'ignore'};
			print OUT "\n";
		}
		my $writing = 0;
		foreach (@infile){
			if (/^BEGIN ORGANISM/){ $writing = 1; next; }
			if ($writing and /^;/){ $writing = 2; }
			next if $writing < 2;
			print OUT $_."\n";
		}
		## show progress ##
		
		if ($iProtein/$nProtein >= $iProgress/60){ 
			print ".";
			$iProgress++;
		}
	}
	
	print " done.\n";
}
print "Execution of filterer completed.\n";
exit 0;

