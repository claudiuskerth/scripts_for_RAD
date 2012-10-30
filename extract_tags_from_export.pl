#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: extract_tags_from_export.pl
#
#        USAGE: ./extract_tags_from_export.pl <file_with_ids> <file_to_extract_from>
#
#  DESCRIPTION: This script is extractimg those loci from an export_sql.pl output
#  				which are found in another file (after filtering, for instance). 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 10/05/12 18:03:39
#     REVISION: ---
#===============================================================================


#
#==========================================================================
#  This program is free software; you can redistribute it and/or modify 
#  it under the terms of the GNU General Public License as published by 
#  the Free Software Foundation; either version 2 of the License, or    
#  (at your option) any later version.                                  
#==========================================================================
#

use strict;
use warnings;

my $usage = "
$0 <file_with_ids> <file_to_extract_from>

This script is intended to extract certain rows from a \".tsv\" sql export file obtained with \"export_sql.pl\"
from a \"stacks\" database. The \"file_with_ids\" should contain the catalog ids in 
its first column and rows with matching catalog ids are extracted from the \"file_to_extract_from\".
The script prints to STDOUT, i. e. pipe into file!
The file with ids should be in the format of the output from stacks export_sql.pl.
The script will extract the catalog ids from this file, so you don't have to prepare
a whitelist file.\n
";

die $usage unless (@ARGV == 2); 
print $usage if $ARGV[0] =~ /-h|-help/;

# get a list of catalog ids which you want to extract from another export_sql.pl output file:
my $file1 = shift;
# skip the header line and cut out the first column
system("tail -n +2 $file1 | cut -f 1 > temp1") == 0 or die $usage;

open(IN1, "<temp1") or die $usage;
my @IDs = <IN1>; # slurp in the whole file
close IN1;
system("rm temp1") == 0 or die $usage;

# open the file which you want to extract from
my $file2 = shift;
open(IN2, "<$file2") or die $usage;

# get the first input line which should contain the header line
my $line = <IN2>;
# split the fields into an array
my @line = split(/\t/, $line);
# get the number of columns for the next split commands
my $col_num = @line;

# print the header line
print join("	", @line);

my %tags;
# store the data in a hash with catalog locus ids as keys
while( ($line = <IN2>) !~ /^$/) {
	@line = split(/\t/, $line, $col_num);
	$tags{$line[0]} = join("	", @line[1..$#line]);
}
close IN2;

# print out only those loci which are present within the first input file 
# (the one containing the "whitelist")
foreach my $ID ( @IDs ) {
	chomp $ID;
	print $ID, "\t", $tags{$ID};
}

exit;
 
