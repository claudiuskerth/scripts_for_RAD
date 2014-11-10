#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: add_RG_tag.pl
#
#        USAGE: ./add_RG_tag.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 09/11/14 20:00:35
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

die "no arguments given\n" unless @ARGV;

my ($file, $RG_given, $RG_present, $RG, @input_files, $outfile);

parse_command_line();

die "no input file given\n" unless @input_files;

if($RG_given){ $RG = $RG_given; }

foreach $file (@input_files){
	if($file =~ /sam$/){ open(IN, $file) or die $!; }
	elsif($file =~ /bam$/){ open(IN, "samtools view -h $file |") or die $!; }
	else{ print "$file does not end in \"sam\" or \"bam\"\n"; next; }

	($outfile = $file) =~ s/(.*)\.(bam|sam)/$1_withRGt.$2/;
	if($outfile =~ /bam$/){ open(OUT, "| samtools view -bhS - > $outfile") or die $!; }
	else{ open(OUT, ">$outfile") or die $!; }


	$RG_present = 0;
	$_ = <IN>;
	while(/^\@/){ 
		print OUT; 
		if(/^\@RG/){ $RG_present = 1; print "RG tag in header. Leaving untouched.\n"; }
		$_ = <IN>; 
	}

	unless($RG_present){
		if($RG_given){ 
			print OUT "\@RG\tID:$RG", "\n";
		}else{
			($RG = $file) =~ s/\.sam|\.bam//;
			$RG =~ s/.*_(par|ery)/$1/;
			print OUT "\@RG\tID:$RG\tSM:$1", "\n";
		}
	}

	$_ = <IN>;
	if(/\tRG:Z:\w+/){ 
		print "Found RG tag in first SAM record. Going to leave RG tags untouched.\n"; 
		print OUT;
		while(<IN>){ print OUT; }
	}else{
		chomp;
		print OUT $_, "\t", "RG:Z:$RG", "\n";
		while(<IN>){
			chomp;
			print OUT $_, "\t", "RG:Z:$RG", "\n";
		}
	}
}

sub parse_command_line {
	while(@ARGV){
		$_ = shift @ARGV;
		if(/^-+RG/i){ $RG_given = shift @ARGV; }
		else{ push @input_files, $_; }
	}	
}
