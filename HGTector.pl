#!/usr/bin/perl

# HGTector (version 0.1.9): Genome-wide detection of horizontal gene transfer based on BLAST hit distribution statistics
# Copyright (C) 2013-2015, Qiyun Zhu, Katharina Dittmar. All rights reserved.
# Licensed under BSD 2-clause license.

use warnings;
use strict;
$| = 1;


## global variables ##

my $i; my $s; my @a;

## welcome information ##

print "
┌---------------------┐
|   HGTector v0.1.9   |
└---------------------┘
";

## check working directory ##

my $scripts = $0;
$scripts =~ s/HGTector\.pl$/scripts/;
my $wkDir = $0;
$wkDir =~ s/HGTector\.pl$/sample/;

if (@ARGV){
	$wkDir = $ARGV[0];
}else{
	print "
Usage:
  1) Select a working directory.
  2) Prepare input protein sequences.
  3) Create a configuration file.
  4) Execute:
     perl HGTector.pl /path/to/working/directory

";
	print "Working directory is not specified. Run a sample task (yes/NO)? ";
	$s = <STDIN>; chomp $s;
	exit 0 if (not $s or $s =~ /^no$/i or $s =~ /^n$/i);
	exit 0 unless ($s =~ /^yes$/i or $s =~ /^y$/i);
	print "Please find the input files under sample/.\n";
	print "Press Enter to proceed, or Ctrl+C to exit:";
	$s = <STDIN>;
}

print "Validating task...\n";
die "Error: Invalid working directory $wkDir.\n" unless -d $wkDir;
print "Working directory: $wkDir.\n";

## read configuration ##

my $interactive = 1;
my $title = "";
my ($trimSeq, $realign, $buildTree);
if (-e "$wkDir/config.txt"){
	open IN, "<$wkDir/config.txt" or die "Error: Invalid configuration file $wkDir/config.txt.\n";
	while (<IN>){
		s/#.*$//; s/\s+$//; s/^\s+//; next unless $_;
		$interactive = 0 if /^interactive=0$/;
		if (/^title=(.+)$/){
			$title = $1;
			print "Job title: $title.\n";
		}
		$trimSeq = $1 if /^trimSeq=(.+)$/;
		$realign = $1 if /^realign=(.+)$/;
		$buildTree = $1 if /^buildTree=(.+)$/;
	}
}else{
	print "Warning: Configuration file $wkDir/config.txt is missing. HGTector will use default settings.\n";
	if ($interactive){
		print "Press Enter to proceed, or Ctrl+C to exit:";
		$s = <STDIN>;
	}
}

## check input protein sets ##

my %sets = ();
print "Validating input files...\n";
die "Error: Input folder $wkDir/input is missing. Please prepare input files before running HGTector.\n" unless -d "$wkDir/input";
opendir (DIR, "$wkDir/input");
@a = readdir(DIR);
closedir DIR;
foreach my $set (@a){
	next if ($set =~ /^\./);
	next unless -s "$wkDir/input/$set";
	open IN, "<$wkDir/input/$set" or next;
	$set =~ s/\.[^\.]+$//;
	$sets{$set} = 0;
	my $intype = "";
	while(<IN>){
		s/\s+$//; next unless $_; next if /^#/;
		unless ($intype){
			if (length($_) >= 12 and substr($_, 0, 12) eq "LOCUS       "){
				$intype = 'GenBank';
				next;
			}elsif (substr($_, 0, 1) eq ">"){
				$intype = 'FASTA';
			}else{
				$intype = 'plain list';
			}
		}
		if ($intype eq 'GenBank'){
			$sets{$set} ++ if (/^\s+\/protein_id=".+"$/);
		}elsif ($intype eq 'FASTA'){
			$sets{$set} ++ if (/^>/);
		}elsif ($intype eq 'plain list'){
			next if /^\d/;
			$sets{$set} ++;
		}
	}
	close IN;
}
die "Error: No input protein sets found in $wkDir/input/\n" unless %sets;
my $nPt = 0;
$nPt += $sets{$_} for keys %sets;
print "Found ".scalar(keys %sets) ." proteins sets ($nPt proteins):\n";
@a = ();
push @a, $_." (".$sets{$_}.")" for (sort keys %sets);
print "  ".join(", ", @a)."\n";

## check BLAST results ##

my $nSet = 0; my $nBla = 0;
if (-d "$wkDir/blast"){
	print "Reading BLAST results from previous runs...\n";
	foreach my $set (sort keys %sets){
		if (-d "$wkDir/blast/$set"){
			opendir (DIR, "$wkDir/blast/$set");
			@a = grep (/\.bla$/, readdir (DIR));
			close DIR;
			if (@a){
				$nSet ++;
				$nBla += scalar (@a);
			}
		}
	}
	print "$nBla BLAST reports for $nSet protein sets are found.\n";
	if (-d "$wkDir/taxonomy" and -s "$wkDir/taxonomy/ranks.db" and -s "$wkDir/taxonomy/taxa.db" and -s "$wkDir/taxonomy/self.info"){
		if ($interactive and $nPt == $nBla){
			print "Has batch BLAST already completed? (yes/NO) ";
			$s = <STDIN>; chomp $s;
			if ($s =~ /^ye?s?$/i){ goto POSTBLAST; }
		}
		print "Validating BLAST results...\n";
		system "$^X $scripts/validator.pl $wkDir";
		$nBla = 0;
		foreach my $set (sort keys %sets){
			if (-d "$wkDir/blast/$set"){
				opendir (DIR, "$wkDir/blast/$set");
				@a = grep (/\.bla$/, readdir (DIR));
				close DIR;
				$nBla += scalar (@a) if (@a);
			}
		}
		goto POSTBLAST if ($nPt == $nBla and not $interactive);
		print "Resuming batch BLAST...\n";
		system "$^X $scripts/blaster.pl $wkDir";
	}else{
		print "Taxonomy information is not found under $wkDir/taxonomy.\n";
		print "The program will run taxonomer.pl to retrieve taxonomic information and then run filterer.pl to filter out invalid or redundant hits.\n";
		if ($interactive){
			print "Proceed? (YES/no) ";
			$s = <STDIN>; chomp $s;
			if ($s =~ /^no?$/i){
				die "Task aborted.\n";
			}
		}
		system "$^X $scripts/taxonomer.pl $wkDir";
		system "$^X $scripts/filterer.pl $wkDir";
	}
	goto POSTBLAST;
}else{
	system "$^X $scripts/blaster.pl $wkDir";
	goto POSTBLAST;
}

POSTBLAST:
print "BLAST results are ready to be analyzed.\n";
if ($trimSeq or $realign or $buildTree){
	if ($interactive){ print "Press Enter to proceed with batch phylogenetic reconstruction, or Ctrl+C to abort:"; $s = <STDIN>; }
	system "$^X $scripts/treer.pl $wkDir";
}
if ($interactive){ print "Press Enter to proceed with prediction, or Ctrl+C to abort:"; $s = <STDIN>; }
system "$^X $scripts/analyzer.pl $wkDir";
if ($interactive){ print "Press Enter to summarize results, or Ctrl+C to abort:"; $s = <STDIN>; }
system "$^X $scripts/summarizer.pl $wkDir";
print "All steps are completed.\n";
exit 0;

