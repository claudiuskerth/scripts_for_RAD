#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: split_by_barcode.pl
#
#        USAGE: ./split_by_barcode.pl  
#
#  DESCRIPTION: takes a file of barcode<tab>individual_name and a fastq file
#  				and assorts fastq records into separate fastq files according
#  				to barcode
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: The script (or rather the "FileHandle" module) prints a "1" to
#         		STDOUT, seemingly for every record.
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 07/09/12 16:25:22
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

use FileHandle;

my $usage = "usage: $0 <barcode_ind_file> <big_fastq_file > /dev/null

Takes a file of barcode<tab>individual_name and a fastq file
and assorts fastq records into separate fastq files according
to barcode. Output files will have the individual name or ID that was
specified in the <barcode_ind_file>.\n";

if ($ARGV[0] =~ /-h|-help/) { die $usage; }
die $usage unless @ARGV == 2;

# read the barcode_ind file into a hash;
# a barcode can only refer to one individual
# but one individual can have several barcodes
open(BC, "<", $ARGV[0]) or die $!;
my %bc_ind;
my ($bc, $ind);
while(<BC>){
	chomp;
	($bc, $ind) = split; # split line on white space
	$bc_ind{$bc} = $ind;
}
close BC;

#foreach my $barcode (sort { $bc_ind{$a} cmp $bc_ind{$b} } keys %bc_ind){
#	print $barcode, "\t", $bc_ind{$barcode}, "\n";
#}

# open the big raw fastq file that should be split
open(IN, "<", $ARGV[1]) or die $!;

# open a file for writing that will take the fastq records that do not match
# any of  the barcodes in the barcode_ind file
open( DIS, ">", "discarded.fq") or die $!;

# initialize
my ($head, $seq, $qh, $qual, $first_5bp, %fh, $diff, $discard);

# while reading in fastq records
while($head = <IN>){
	$seq = <IN>;
	$qh = <IN>;
	$qual = <IN>;

	$discard = 1;

	# get the first 5bp of the read
	$first_5bp = substr($seq, 0, 5);

	foreach my $barcode (keys %bc_ind){
		# foreach barcode, get the number of mismatches
		$diff = ($barcode ^ $first_5bp) =~ tr[\001-\255][];
		# if the number of mismatches is 0 or 1
		if ($diff <= 1){
			# then the barcode is in the provided list (barcode_ind)
			# and the record won't be discarded
			$discard = 0;
			# get the individual name or id for the barcode
			$ind = $bc_ind{$barcode};
			# if no output file has been opened yet for this individual,
			# open one with the individual name as the file name
			if( not defined $fh{$ind} ){
				$fh{$ind} = FileHandle->new;
				$fh{$ind}->open(">$ind.fq") or die $!;
			}
			# print to the file of the individual to which this barcode belongs
			print $fh{$ind}->print($head, $seq, $qh, $qual);
		}
	}
	# if no matching barcode could be found, dump this record to a separate file
	if ($discard) {print DIS $head, $seq, $qh, $qual;}
}
# close all open file handles
foreach my $ind ( keys %fh ){
	$fh{$ind}->close
}

close DIS;
close IN;

autoflush STDOUT 1;
