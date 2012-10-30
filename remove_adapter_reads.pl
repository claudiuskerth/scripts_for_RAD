#!/usr/bin/env perl
#
#===============================================================================
#
#         FILE: remove_adapter_reads.pl
#
#        USAGE: ./remove_adapter_reads.pl  
#
#  DESCRIPTION: this script takes a fastq file from standard input and outputs
#  				a file to standard output (i. e. pipe into file) 
#  				that should be largely free of reads from adapter dimers
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 24/09/12 17:20:21
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

my $usage = "
this script takes a fastq file from standard input and outputs
a file to standard output (i. e. pipe into file) 
that should be largely free of reads from adapter dimers
\n";

die $usage if( defined($ARGV[0]) && $ARGV[0] =~ /^-h$/ );

my ($head, $seq, $qh, $qual); 

my $adapter_seq1 = "^CGTAT"; # don't use the barcode "CGTAT"
my $adapter_seq2 = "AGATCGGAAG"; # first 10bp of illumina adapters

# end of P2 illumina adapter after selective PCR, 
# this sequence is usually followed by a more or less long Adenine homopolymer:
my $adapter_seq3 = "CGTATGCCGTCTTCTGCTTG"; 

while(<>){
	$head = $_;
	$seq = <>;
	$qh = <>;
	$qual = <>;
	print $head, $seq, $qh, $qual unless ($seq =~ /$adapter_seq1|$adapter_seq2|$adapter_seq3/o); 
}
