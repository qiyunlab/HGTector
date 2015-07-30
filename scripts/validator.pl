#!/usr/bin/perl

use warnings;
use strict;
$| = 1;

print "

This script checks the completeness of batch BLAST results. By default, it deletes
incomplete BLAST reports (*.bla) or invalid hits from each report.

Usage:
  perl validator.pl <working directory>

Note:
  Add parameter \"-keep\" to display a report without deleting invalid results.

Output:
  screen output

" and exit unless @ARGV;

print "
-> Validator: Check completeness of batch BLAST results. <-
";

my $wkDir = $ARGV[0];
my $keep = 0; $keep = 1 if ($#ARGV and $ARGV[1] eq "-keep"); # don't delete invalid reports or hits

opendir (DIR, "$wkDir/blast");
my @sets = grep {-d "$wkDir/blast/$_" and not /^\.{1,2}$/} readdir(DIR);
close DIR;
die "No BLAST results detected.\n" unless @sets;

foreach my $set (@sets){
	print "Checking protein set $set...\n";
	opendir DIR, "$wkDir/blast/$set";
	my @blas = grep (s/\.bla$//, readdir(DIR));
	close DIR;
	my $allHit = 0;
	my $badBla = 0;
	my $badHit = 0;
	foreach my $bla (@blas){
		my @badHits = ();
		my @out = ();
		open IN, "<$wkDir/blast/$set/$bla.bla";
		while (<IN>){ s/\s+$//; push (@out, $_); }
		close IN;
		my $k = 1; # status of report
		foreach (@out){
			next unless $_;
			if ($k == 1 and (/^#NEXUS/)){ $k = 2; next; }
			if ($k == 2 and /^BEGIN QUERY;/){ $k = 3; next; }
			# if ($k == 3 and /^\tGI=\d+;$/){ $k = 4; next; }
			# if ($k == 4 and /^\tAccession=.+;$/){ $k = 5; next; }
			# if ($k == 5 and /^\tLength=\d+;$/){ $k = 6; next; }
			# if ($k == 6 and /^\tProduct=.+;$/){ $k = 7; next; }
			# if ($k == 7 and /^\tOrganism=.+;$/){ $k = 8; next; }
			if ($k == 3 and /^END;$/){ $k = 4; next; }
			if ($k == 4 and /^BEGIN ORGANISM;$/){ $k = 5; next; }
			if ($k == 5 and /^END;$/){ $k = 6; last; }
			if ($k == 5){
				next if /^\[/;
				next if /^;/;
				my @b = split (/\t/);
				if (@b and $#b < 5){
					$_ = "[deleted]";
					push (@badHits, $b[0]);
				}
				$allHit ++;
				next;
			}
		}
		if ($k < 6){ # bad report
			$badBla ++;
			print "Incomplete BLAST report for $bla";
			unless ($keep){ unlink "$wkDir/blast/$set/$bla.bla"; print ". Deleted"; }
			print ".\n";
		}else{
			if (@badHits){
				$badHit += scalar (@badHits);
				print "Invaid hits for $bla:";
				print " $_" for (@badHits);
				unless ($keep){
					open OUT, ">$wkDir/blast/$set/$bla.bla";
					foreach (@out){
						next if /^\[deleted\]/;
						print OUT "$_\n";
					}
					close OUT;
					print ". Deleted";
				}
				print ".\n";
			}
		}
	}
	print "Protein set $set checked. $badBla of ".(scalar @blas)." BLAST reports are incomplete. $badHit of $allHit hits are invalid.\n";
}
print "done.\n";


