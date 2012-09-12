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

my $usage = "usage: $0 [-s] -bc <barcode_ind_file> -i <big_fastq_file > /dev/null

Takes a file of barcode<tab>individual_name and a fastq file
and assorts fastq records into separate fastq files according
to barcode. Barcodes need to have equal length. By default, the script allows
up to 1bp mismatch to a barcode in the <barcode_ind> file, so 3bp mismatch barcodes
are recommended to avoid ambiguous assignments. Output files will 
have the individual name or ID that was specified in the <barcode_ind_file>.
This script has the feature that it can be used when you have used more than
one barcode for an individual in order to diversify the first few base pairs
in the read for the illumina sequencer. Then, all records with barcodes matching
a certain individual will be printed into one file. The barcodes can be stripped off
in the process. If you have split by barcode with another programme that creates a
separate file for each barcode, but multiple barcodes were used to identify one
individual, then use my script \"rename_sample_files.pl\".

options:

-s turns on stripping of barcodes from reads
-bc barcode_ind_file
-i fastq file to split by barcode
\n";

my ($strip_barcode, $barcode_ind_file, $fastq_in);

# parse command line
while (@ARGV){
	$_ = shift;
	if (/-s/) { $strip_barcode = 1; }
	if (/-bc/) { $barcode_ind_file = shift; }
	if (/-i/) { $fastq_in = shift; }
	if (/-h|-help/) { die $usage; }
}

die $usage unless ( defined($barcode_ind_file) && defined($fastq_in) );

# read the barcode_ind file into a hash;
# a barcode can only refer to one individual
# but one individual can have several barcodes
open(BC, "<", $barcode_ind_file) or die $!;
my %bc_ind;
my ($bc, $ind, $bc_length);
while(<BC>){
	chomp;
	($bc, $ind) = split; # split line on white space
	$bc_ind{$bc} = $ind;
}
$bc_length = length($bc);
close BC;

#foreach my $barcode (sort { $bc_ind{$a} cmp $bc_ind{$b} } keys %bc_ind){
#	print $barcode, "\t", $bc_ind{$barcode}, "\n";
#}

# open the big raw fastq file that should be split
open(IN, "<", $fastq_in) or die $!;

# open a file for writing that will take the fastq records that do not match
# any of  the barcodes in the barcode_ind file
open( DIS, ">", "discarded.fq") or die $!;

# initialize
my ($head, $seq, $qh, $qual, $first_5bp, %fh, $diff, $discard, $bc_from_read);

# while reading in fastq records
while($head = <IN>){
	$seq = <IN>;
	$qh = <IN>;
	$qual = <IN>;

	$discard = 1;

	# get the barcode from the read
	$bc_from_read = substr($seq, 0, $bc_length);

	foreach my $barcode (keys %bc_ind){
		# foreach barcode, get the number of mismatches
		$diff = ($barcode ^ $bc_from_read) =~ tr[\001-\255][];
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
			# strip barcode sequence of read if specified on the command line
			if ($strip_barcode){
				$seq = substr($seq, $bc_length);
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
