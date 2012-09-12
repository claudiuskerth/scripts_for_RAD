#!/usr/bin/env perl
use strict; use warnings;

my $usage = "\n$0 [options] <barcode_ind file> <file names to rename> 

This script is designed to rename output files from stacks' \'process_radtags\'script, 
which contain the barcode sequence in their file names. This script replaces the
barcode sequence in the filenames with a name of your liking.
<file names to rename> can include a wildcard expression which needs to be in double quotes,
for instance \"*.fq_1\".
<barcode_ind file> should be a tab delimited file in the format of:
			<barcode>	<individual name>
Example:
			AGCTA	30-1
			ATCGT	30-2
			etc.
If more than one barcode matches an <individual name> in the <barcode_ind file>, 
then all sequence records for those barcodes will be concatenated into one file.
This can be useful if you have used more than one barcode per individual in order to diversify
the beginning of the reads for the illumina sequencer.

-n no action, just print what would be done 
\n";

die $usage unless @ARGV >= 2;

my @files = glob( pop(@ARGV) ) or die $!;

#foreach my $file (@files) { print $file, "\n"; }

open(IN, "<", pop(@ARGV) ) or die $!;

my ($barcode, $ind, %barcode_ind, $no_action);

if ( defined($ARGV[0]) && $ARGV[0] =~ /-n/) { $no_action = "true"; }

while(<IN>){
	chomp;
	($barcode, $ind) = split /\t/, $_;
	next unless ( defined($barcode) && defined($ind) );
	$barcode_ind{$barcode} = $ind;
}
close IN;
#foreach my $barcode (sort keys %barcode_ind) {
#print $barcode, "\t", $ind, "\n";
#}

File: 
# foreach file to be renamed
foreach my $filename (@files) {
	# look for the individual name this barcode stands for 
	foreach my $barcode ( sort keys %barcode_ind ) {
		if ($filename =~ /$barcode/) {
			# get a new file name by replacing the barcode sequence with the
			# individual name
			(my $new_filename = $filename) =~ s/$barcode/$barcode_ind{$barcode}/;
			# if the new file name already exists because you have used more than
			# one barcode per individual in order to diversify the start of the reads
			if (-e $new_filename) {
				# add the contents of this file to the end of the already existing file
				system("cat $filename >> $new_filename") == 0 or die "cannot concatenate\n";
				system("rm $filename") == 0 or die $!;
			} else {
				if ($no_action) { print "$filename will be renamed as $new_filename\n"; }
				# rename the file with the new file name
				else{ system("mv $filename $new_filename") == 0 or die "cannot rename\n"; }
			}
			next File;
		}
	}
}
