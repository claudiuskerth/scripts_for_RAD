#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: SAM_flag_dist.pl
#
#        USAGE: ./SAM_flag_dist.pl  
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
#      CREATED: 31/10/16 17:23:28
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

my $usage = "
$0 <BAM/SAM file name>
\n";

# check that samtools is in PATH
my $samtools = `which samtools`;
die "need samtools in PATH" unless $samtools =~ /samtools/;

# die if no BAM or SAM file name provided
die $usage unless $ARGV[0];

# open stream of SAM flags
open( my $IN, "samtools view $ARGV[0] | cut -f 2 |") or die $!;

my (%flag, $k, $v);

# tally SAM flags
while(<$IN>){
	chomp;
	$flag{$_}++;
}

# pipe output into sort
open( my $OUT, "| sort -gr -k 2,2") or die $!;

# print out each flag, count and flag description
while( ($k, $v) = each %flag ){
	printf $OUT "%3d\t%7d\t%s", $k, $v, `samtools flags $k | cut -f3`;
}
