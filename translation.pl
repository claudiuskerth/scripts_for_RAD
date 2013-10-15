#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: translation.pl
#
#        USAGE: ./translation.pl  
#
#  DESCRIPTION: This script takes a (multi-) fasta file of DNA sequences as 
#  				input and prints info about the longest open reading frame in each.
#  
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), claudiuskerth[at]gmail.com
# ORGANIZATION: Sheffield University
#      VERSION: 1.0
#      CREATED: 15/10/2013 17:24:54
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

# check whether required module is installed
eval { require Bio::SeqIO };
if ($@) { 
	die "$0 requires the Bioperl module Bio::SeqIO. 
Please install this package and add to your Perl library path.\n";
}else{ Bio::SeqIO->import() };

# usage statement
my $usage = "
This script takes a (multi-) fasta file of DNA sequences as 
input and prints info about the longest open reading frame in each.
\n$0 <fasta file to translate>\n\n";

die $usage unless @ARGV;

my $infile = $ARGV[0];

# create object for file reading
my $seqIO = Bio::SeqIO->new(-file => $infile,
							-format => "fasta"
							);
							
# initialise variables
my ($rf_1,  $rf_2,  $rf_3,  $rf_4,  $rf_5,  $rf_6,  @rf_1,  @rf_2,  @rf_3,  @rf_4,  @rf_5,  @rf_6,  %hash);
my $longest_ORF;

# foreach fasta sequence:
while(my $seq_obj = $seqIO->next_seq){
	print $seq_obj->display_id, "\n";
	print $seq_obj->desc, "\n";
	# get translation in reading frame 1
	$rf_1 = $seq_obj->translate(-frame => 0)->seq; 
	# split the translation at stop codons
	@rf_1 = split('\*', $rf_1); 
	# store the longest putative ORF of this frame in a hash
	$hash{ longest_element(@rf_1) } = "reading frame 1"; 
	#foreach(@rf_1){ print length($_), "\n";}
	# get translation in reading frame 2
	$rf_2 = $seq_obj->translate(-frame => 1)->seq; 
	@rf_2 = split('\*', $rf_2);
	$hash{ longest_element(@rf_2) } = "reading frame 2";
	$rf_3 = $seq_obj->translate(-frame => 2)->seq;
	@rf_3 = split('\*', $rf_3);
	$hash{ longest_element(@rf_3) } = "reading frame 3";
	$rf_4 = $seq_obj->revcom->translate(-frame => 0)->seq;
	@rf_4 = split('\*', $rf_4);
	$hash{ longest_element(@rf_4) } = "reading frame 4";
	$rf_5 = $seq_obj->revcom->translate(-frame => 1)->seq;
	@rf_5 = split('\*', $rf_5);
	$hash{ longest_element(@rf_5) } = "reading frame 5";
	$rf_6 = $seq_obj->revcom->translate(-frame => 2)->seq;
	@rf_6 = split('\*', $rf_6);
	$hash{ longest_element(@rf_6) } = "reading frame 6";

	# get the longest ORF of all reading frames as protein sequence
	$longest_ORF = longest_element( keys(%hash) ); 
	# print the reading frame which produced the longest ORF
	print $hash{ $longest_ORF }, "\n"; 
	# print out protein sequence of longest ORF along 
	# with accession and description form the fasta header
	print $longest_ORF, "\t", $seq_obj->display_id, "\t", $seq_obj->desc, "\n"; 	
	print length( $longest_ORF ) ;
	print "\n\n";
}
#####################################
# SUBROUTINE
# name: longest_element
# receives: array of wutative ORF's
# returns: longest ORF in array
# from http://stackoverflow.com/questions/4182010/the-fastest-way-execution-time-to-find-the-longest-element-in-an-list
#####################################
sub longest_element{
	my $max = -1;
	my $max_ref;
	for(@_){
		if(length($_) > $max){
			$max = length($_);
			$max_ref = \$_;
		}
	}
	return $$max_ref;
}
#####################################
