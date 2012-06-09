#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: remove_dodgy_positions.pl
#
#        USAGE: ./remove_dodgy_positions.pl  
#
#  DESCRIPTION: This script is designed to remove certain positions from all reads
#  				in the input files. Currently this file will remove the 96th and 196th
#  				position.
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 01/06/12 12:51:03
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

my $usage = "
$0 [-h] [input_files]

This script is designed to remove certain positions from all reads
in the input files. Currently this file will remove the 96th and 196th
position. The input files can be provided as a glob like \"*.fq\".
Input files need to end with \".fq\".

\n";

foreach(@ARGV){
	if ($_ =~ /-h/){ die $usage; }
}

die $usage unless @ARGV>0;

# initialize
my ($header, $seq, $qual_head, $qual, $seq_a, $seq_b, $qual_a, $qual_b);

# iterate over input file (can be a glob expanded by bash, like *.fq)
foreach my $file (@ARGV){

	# open input file or reading
	open(IN, "<$file") or die $!;
	
	# make output file name
	$file =~ s/\.fq$/_spliced.fq/;

	# open output file
	open(OUT, ">$file") or die $!;
	print STDERR "writing to file ... $file\n";
	
	# iterate over input file
	while(<IN>){
		# read in one record
		$header = $_;
		$seq = <IN>;
		$qual_head = <IN>;
		$qual = <IN>;
	
		# get the first 95 bases
		$seq_a = substr($seq, 0, 95);
		# get the 97th to 195th base of the merged reads
		$seq_b = substr($seq, 96, 99);
		# see file /data/claudius/DD-RAD/second_run/export/20052012/SNPs_across_tag_pos
	
		$qual_a = substr($qual, 0, 95);
		$qual_b = substr($qual, 96, 99);

		# print the spliced record to output file
		print OUT $header,
				  $seq_a . $seq_b . "\n",
				  $qual_head,
				  $qual_a . $qual_b . "\n";
	}
	close IN;
	close OUT;
}
exit;


