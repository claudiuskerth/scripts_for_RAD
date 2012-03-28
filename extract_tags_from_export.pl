#!/usr/bin/env perl
use strict; use warnings;

my $usage = "\n$0 <file_with_ids> <file_to_extract_from>\n\n" . 
"prints to STDOUT, i. e. pipe into file!\n" .
"The file with ids should be in the format of the output from stacks export_sql.pl.\n" .
"This script will extract the catalog ids from this file, so you don't have to prepare\n" .
"a whitelist file.\n\n";

die $usage unless (@ARGV == 2); 

# get a list of catalog ids which you want to extract from another export_sql.pl output file:
my $file1 = shift;
# skip the header line and cut out the first column
system("tail -n +2 $file1 | cut -f 1 > temp1") == 0 or die $usage;

open(IN1, "<temp1") or die $usage;
my @IDs = <IN1>;
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
# store the data in hash with catalog locus ids as keys
while($line = <IN2>) {
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
 
