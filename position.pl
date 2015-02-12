#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: position.pl
#
#        USAGE: ./position.pl  
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
#      CREATED: 04/02/15 15:12:41
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

my $file;
my %SbfI;

foreach my $filename (@ARGV){
	# open file with filename, unzip if necessary
	if($filename =~ /gz$/){
		open( $file, "zcat $filename |") or die $!;
	}else{
		open( $file, $filename) or die $!;
	}
	my %pos;
	# read in the file line by line
	while(<$file>){
		# skip lines that do not contain sequence
		next if ($.-2)%4;
		# if the sequence contains an SbfI site, then
		# store the 1-based position in the read as a value
		# in a hash with the read sequence as key. That way
		# identical sequences are recorded only once.
		if(/CCTGCAGG/){ $pos{$_} = $-[0]+1 }
	}
	$SbfI{$filename} = join(",", values %pos);
}

foreach my $filename (keys %SbfI){
	print $filename, ",", $SbfI{$filename}, "\n";
}

