#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: decode_SAMflag.pl
#
#        USAGE: ./decode_SAMflag.pl  
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
#      CREATED: 24/10/16 21:57:19
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

my $usage = "
$0 <SAM flag [INT]>
\n";

die $usage unless @ARGV;

my $flag = $ARGV[0];

die $usage if $flag =~ /\D/;


if($flag & 1){ print "read from paired-end sequencing\n"; }
else{ print "read from single-end sequencing\n"; }
if($flag & 2){ print "read is mapped in proper pair :-)\n"; }
else{ print "read not mapped in proper pair :-(\n"; }
if($flag & 4){ print "read is unmapped\n"; }
else{ print "read is mapped :-)\n"; }
if($flag & 8){ print "mate is unmapped :-(\n"; }
else{ print "mate is mapped :-)\n"; }
if($flag & 16){ print "read sequence printed as revcomp of input\n"; }
# Unmapped reads can have strand. If true this flag only means the sequence
# is printed as reverse complement in the SAM file.
else{ print "read sequence printed as input\n"; }
if($flag & 32){ print "mate sequence printed as revcomp of input\n"; }
else{ print "mate sequence printed as input\n"; }
if($flag & 64){ print "read is the first read in a pair\n"; }
if($flag & 128){ print "read is the second read in a pair\n"; }
if($flag & 256){ print "read is not primary (a read having split hits may have multiple primary alignment records)\n"; }
if($flag & 512){ print "read fails platform/vendor quality checks\n"; }
if($flag & 1024){ print "read is PCR or optical duplicate\n"; }
if($flag & 2048){ print "supplementary alignment\n"; }
