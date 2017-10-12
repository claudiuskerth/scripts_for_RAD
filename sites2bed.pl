#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: sites2bed.pl
#
#        USAGE: ./sites2bed.pl  
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
#      CREATED: 22/11/16 15:17:54
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

$_ = <>; # read in first line of sites file
my ($previous_contig, $previous_pos) = split; # split on white space and record contig name and first position
my $begin = $previous_pos - 1; # record beginning of first interval, 0-based coordinates

my $current_contig = "";
my $current_pos = 0;

while(<>){
	($current_contig, $current_pos) = split;	
	# an interval needs to be a sequence of contiguous sites. If either the contig name changes
	# or the new position is not 1 above the previous, then the current interval ends with the previous
	# conbtig and previous position and a new interval shall be started:
	if($current_contig ne $previous_contig or $current_pos != ($previous_pos + 1)){
		print $previous_contig, "\t", $begin, "\t", $previous_pos, "\n"; # print out interval in BED format
		$begin = $current_pos - 1; # initialise beginning of new interval
		$previous_contig = $current_contig;
	}
	$previous_pos = $current_pos;
}
# if end of file is reached, print out the last interval
print $previous_contig, "\t", $begin, "\t", $previous_pos, "\n";
