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
\n$0 \n
	-collect  --  print sequences with identical ORF's to the same file
	-in      <fasta file to translate>
\n";

die $usage unless @ARGV;

my ($infile, $collect_identical_ORF);

parse_command_line();


# create object for file reading
my $file = Bio::SeqIO->new(-file => $infile,
							-format => "fasta"
							);
							
# initialise variables
my ($rf_1,  $rf_2,  $rf_3,  $rf_4,  $rf_5,  $rf_6,  @rf_1,  @rf_2,  @rf_3,  @rf_4,  @rf_5,  @rf_6);
my $longest_ORF;
my %protein_seq;

# foreach fasta sequence:
while(my $fasta = $file->next_seq){
	my %ORF;
	print $fasta->length, "\n";
	# get translation in reading frame 1
	$rf_1 = $fasta->translate(-frame => 0)->seq; 
	# split the translation at stop codons
	@rf_1 = split('\*', $rf_1); 
	# store the longest putative ORF of this frame in a hash
	$ORF{ longest_element(@rf_1) } = "reading frame 1"; 
	#foreach(@rf_1){ print length($_), "\n";}
	# get translation in reading frame 2
	$rf_2 = $fasta->translate(-frame => 1)->seq; 
	@rf_2 = split('\*', $rf_2);
	$ORF{ longest_element(@rf_2) } = "reading frame 2";
	$rf_3 = $fasta->translate(-frame => 2)->seq;
	@rf_3 = split('\*', $rf_3);
	$ORF{ longest_element(@rf_3) } = "reading frame 3";
	$rf_4 = $fasta->revcom->translate(-frame => 0)->seq;
	@rf_4 = split('\*', $rf_4);
	$ORF{ longest_element(@rf_4) } = "reading frame 4";
	$rf_5 = $fasta->revcom->translate(-frame => 1)->seq;
	@rf_5 = split('\*', $rf_5);
	$ORF{ longest_element(@rf_5) } = "reading frame 5";
	$rf_6 = $fasta->revcom->translate(-frame => 2)->seq;
	@rf_6 = split('\*', $rf_6);
	$ORF{ longest_element(@rf_6) } = "reading frame 6";
	
	# get the longest ORF of all reading frames as protein sequence
	$longest_ORF = longest_element( keys(%ORF) ); 
	# print report
	print $fasta->display_id, "\t", $ORF{ $longest_ORF },"\n";
	print $fasta->desc, "\n";
	# print the reading frame which produced the longest ORF
	print $ORF{ $longest_ORF }, "\n"; 
	# print out protein sequence of longest ORF along 
	# with accession and description form the fasta header
	print $longest_ORF, "\t",  $ORF{ $longest_ORF }, "\n"; 	
	print length( $longest_ORF ) ;
	print "\n\n";

	if($collect_identical_ORF){
		# store the sequence object in an anonymous array, which is the
		# value to the longest ORF protein sequence of the DNA sequence as key
		push @{ $protein_seq{$longest_ORF} },  $fasta;	
	}
}

if($collect_identical_ORF){
	# print fasta records with identical ORF protein sequences 
	# together into one file
	my $i = 1;
	foreach my $longest_ORF ( keys %protein_seq ){
		# create object for file writing
		my $outfile= Bio::SeqIO->new(-file => ">$i.fasta",
									-format => "fasta"
									);
		foreach my $fasta ( @{ $protein_seq{$longest_ORF} } ){
			$outfile->write_seq($fasta);
		}
		$i++;
	}
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

sub parse_command_line{
	while(@ARGV){
		$_ = shift @ARGV;
		if(/^-collect$/){$collect_identical_ORF = 1}
		if(/^-in$/){$infile = shift @ARGV}
		die "\nno input file given\n" unless $infile;
	}
}
