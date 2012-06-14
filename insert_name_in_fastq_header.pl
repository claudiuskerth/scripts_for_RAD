#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: insert_name_in_fastq_header.pl
#
#        USAGE: ./insert_name_in_fastq_header.pl  
#
#  DESCRIPTION: This script inserts the sample name which is part of the input
#  				file name into the fastq header right after the barcode sequence
#  				(inserted by 'process_radtags').
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 12/06/12 18:44:01
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

# foreach input file given on the command line (can be glob expanded by the shell):
foreach(@ARGV){
	# open the file for reading
	open(IN, "<", $_) or die $!;

	# remove everything from the second underscore toward the end of the file name
	# that should be the sample name
	s/^(.+?_.+?)_.+$/$1/;
	my $name = $_;

	# open output file
	open(OUT, ">", "$name" . "_new_fq_header.fq") or die $!;

	# read in the file line by line
	while(<IN>){
		# if line is fastq header from 'process_radtags'
		if ( /^(@[ATGC]{5}_)(.*$)/ ){
			# insert sample name
			print OUT $1, "_", $name, "__$2\n";
		}
		else { print OUT; }
	}
	close IN;
	close OUT;
}
