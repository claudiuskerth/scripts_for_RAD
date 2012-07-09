#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: reformat_qual_scores.pl
#
#        USAGE: ./reformat_qual_scores.pl  
#
#  DESCRIPTION: This script is designed to convert the Phred+33 quality score
#  				characters in a .fastq file into CAP3/Phrap compatible numbers.
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 13/06/12 17:29:13
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

my $usage = "
$0 <input files>

This script is designed to convert the Phred+33 quality score
characters in a .fastq file into CAP3/Phrap compatible numbers.
It produces two output files for each input file, one sequence file
and one quality file (both in fasta format).
\n";


die $usage unless (@ARGV>0);

# foreach input file given on the command line (can be a glob expanded by the shell)
foreach(@ARGV){
	# open the file for reading
	open(IN, "<", $_) or die $!;
	# replace file ending
	s/\.fq$/.fasta/;
	# open output file for fasta sequence
    open(OUT_SEQ, ">", $_) or die $!;
#	print $_, "\n";
#	# replace file ending
#	s|\.fasta|\.qual|;
	# open output file for quality scores
	open(OUT_QUAL, ">", "$_.qual") or die $!;
#	print $_, "\n";
	# read in input file line by line
	while(<IN>){
		my $header = $_;
		# make fasta header from fastq header 
		$header =~ s/^@/>/;
		my $seq = <IN>;
		# print header and sequence to output in fasta format
		print OUT_SEQ $header, $seq;
		my $qh = <IN>;
		my $qual = <IN>;
		chomp $qual;
		# replace quality score ASCII characters with their numerical meaning (Phred+33 format)
		$qual =~ s/(.)/sprintf("%s ", ord($&)-33)/eg; # behold the power of perl's substitution command !!!
		print OUT_QUAL $header, $qual, "\n";
	}
	close IN;
	close OUT_SEQ;
	close OUT_QUAL
}

