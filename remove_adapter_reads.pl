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
./remove_adapter_reads.pl [options] <inputfile>

this script takes a fastq file from standard input and outputs
a file to standard output (i. e. pipe into file) 
that should be largely free of reads from adapter dimers

Command line option:

	-fuzzy [number mismatches] : enable fuzzy matching and specify the degree of fuzziness
	-ac [1|2] : choose between search sequences:
 	   1: AGATCGGAAG (default)
	   2: CGTATGCCGTCTTCTGCTTG
	-adapter : provide custom search sequence (i. e. adapter seq)
\n";

my $fuzzy = 0;

my $fuzziness = 0;

my $adapter_seq = "";

my $adapter_seq1 = "AGATCGGAAG"; # first 10bp of illumina adapters

# end of P2 illumina adapter after selective PCR, 
# this sequence is usually followed by a more or less long Adenine homopolymer:
my $adapter_seq2 = "CGTATGCCGTCTTCTGCTTG"; 

my $adapter_choice = 1;

parse_command_line(); 

if($adapter_choice == 1){ $adapter_seq = $adapter_seq1; }
elsif($adapter_choice == 2){ $adapter_seq = $adapter_seq2; }
else{ die "Please pick adapter 1 or 2 pr provide your own adapter sequence to the \"-ac\" switch.\n" }

my $al = length($adapter_seq);
if($fuzziness >= $al){ die "Fuzziness, i. e. allowed number of mismatches, is greater than search sequence!\n" }

# check
#print "$fuzzy, $fuzziness, $adapter_seq, $adapter_choice\n";
#exit;

my ($head, $seq, $qh, $qual); 

#my $adapter_seq1 = "^CGTAT"; # don't use the barcode "CGTAT"


if($fuzzy){
	my $window = "";
	my $sl = 0;
	my $mismatch = 0;
RECORD: while(<>){
		$head = $_;
		$seq = <>;
		$qh = <>;
		$qual = <>;
		$sl = length($seq);
		for(my $i=0; $i<$sl-$al; $i++){
			$window = substr($seq, $i, $al);	
			$mismatch = ($window ^ $adapter_seq) =~ tr/\001-\255//;
			if ($mismatch <= $fuzziness){
				print "Found one!\n";
				next RECORD; 
			}
		}
		print $head, $seq, $qh, $qual; 
	}
}else{
	while(<>){
		$head = $_;
		$seq = <>;
		$qh = <>;
		$qual = <>;
		print $head, $seq, $qh, $qual unless ($seq =~ /$adapter_seq1|$adapter_seq2/o); 

	}
}

sub parse_command_line{
	die $usage unless @ARGV;
	while(@ARGV>1){
		$_ = shift(@ARGV);
		if(/^-fuzzy$/){ $fuzzy = 1; $fuzziness = shift @ARGV; }
		if($fuzziness =~ /[^0-9]/){ die "Please specify fuzzyness.\n" }
		if(/^-adapter$/){ $adapter_seq = shift @ARGV; }
		if(/^-ac$/){ $adapter_choice = shift @ARGV; }
		if(/^-h$/){ die $usage; }
	}
}
