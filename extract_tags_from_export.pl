#!/usr/bin/env perl
use strict; use warnings;

my $usage = "./extract_tags_from_export.pl <file_with_ids> <file_to_extract_from>
prints to STDOUT, i. e. pipe into file!
"

unless (@ARGV == 2) die $usage;
# get a list of catalog ids which you want to extract from another export_sql.pl output file:
my $file1 = shift;

system("tail -n +2 $file1 | cut -f 1 > temp1") == 0 or die $usage;

open(IN1, "<temp1") or die $usage;
my @IDs = <IN1>;
close IN1;
system("rm temp1") == 0 or die $usage;

my $file2 = shift;

# open the file from which you want to extract from
open(IN2, "<$file2") or die $usage;

# get the first input line which should contain the header line
my $line = <IN2>;
# split the fields into an array
my @line = split(/\t/, $line);
# get the number of columns for the next split commands
my $col_num = @line;

print join("	", @line);

my %tags;

while($line = <IN2>) {
	@line = split(/\t/, $line, $col_num);
	$tags{$line[0]} = join("	", @line[1..$#line]);
}
close IN2;

foreach my $ID ( @IDs ) {
	chomp $ID;
	print $ID, "\t", $tags{$ID};
}


 
