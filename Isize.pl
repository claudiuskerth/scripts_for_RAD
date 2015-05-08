#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: Isize.pl
#
#        USAGE: ./Isize.pl  
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
#      CREATED: 07/05/15 09:36:57
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

my (@F, $c, @M);

while(<>){
	next if /^\@/;
	@F = split;
	if($F[1] & 16){
		@M = $F[5] =~ m/(\d+)M/g;
		$c = 0;
		for(@M){ $c += $_ };
		print join("\t", @F[0..8]), "\t";
		print $F[3] + $c - $F[7], "\t";
		print join("\t", @F[9..$#F]), "\n";
	}else{
		print join("\t", @F[0..8]), "\t\t";
		print join("\t", @F[9..$#F]), "\n";
	}
}
