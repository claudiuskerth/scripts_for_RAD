#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: sort_by_column.pl
#
#        USAGE: ./sort_by_column.pl  
#
#  DESCRIPTION: sorts ascending by specified column
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 20/05/12 22:42:49
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

my $column=1;
die "no argument\n" unless @ARGV;

$_ = shift @ARGV;
if ($_ =~ /-n/ ) { $column = shift @ARGV; }
else { die "please specify the column according to which the file should be sorted.\n"; } 
#
#my $infile = shift @ARGV;
#chomp $infile;

#open(IN, "<STDIN") or die $!;

my %hash;

while(<>){
	last if $_ =~ /^$/;
	chomp;
	my @line = split(/\t/, $_);
	$hash{$line[$column-1]} = $_;
}

foreach my $field ( sort {$a <=> $b} keys %hash ) {
	print $hash{$field}, "\n";
}


