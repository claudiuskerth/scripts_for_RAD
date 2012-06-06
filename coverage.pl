#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: coverage.pl
#
#        USAGE: ./coverage.pl  <input_file_with_allele_depths>
#
#  DESCRIPTION: Calculates average coverage per allele for each locus in an export
#  				output file from export_sql.pl containing allele depths (instead of
#  				genotypes). 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 10/05/12 18:19:07
#     REVISION: ---
#===============================================================================

use strict;
use warnings;


my $usage = "
$0 <input_file_with_allele_depths>

Prints to STDOUT, i. e. pipe into file!
This script takes an sql export output file containing allele depths created with \"export_sql.pl\" from a \"stacks\"
database and calculates average coverage per allele for each locus over the whole sample. This is reported in an extra
column that is inserted in the output table.
\n";

die $usage unless @ARGV == 1;
die $usage if $ARGV[0] =~ /-h|-help/;

# get the file name to read
my $infile = shift or die;

# open file to read
open(IN, "<$infile") or die;

# read in the first line
my $line = <IN>;
my @line = split(/\t/, $line);
my $num_col = @line;

print STDERR "please type in the name of the leftmost sample in the table (at least its unique part):\n";
my $first_sample_name = <STDIN>;
chomp $first_sample_name;

print STDERR "please type in the name of the rightmost sample in the table (at least its unique part):\n";
my $last_sample_name = <STDIN>;
chomp $last_sample_name;

# get the index for the column of the leftmost sample in the table 
my ($first_sample_col) = grep { $line[$_] =~ /$first_sample_name/  } 0..$#line;

# get the index for the column of the rightmost sample in the table 
my ($last_sample_col) = grep { $line[$_] =~ /$last_sample_name/  } 0..$#line;

# print headers and add two new columns
print join("	", @line[0..$first_sample_col-1], 
		"geno_calls", 
		"average_coverage/allele", 
		"total_coverage/locus",
		 @line[$first_sample_col..$#line]);

open(ERR, ">coverage.err") or die $!;

my ($sum_cov, $mean_cov, $geno_calls, $allele_count);

# while reading in loci
while ($line = <IN>) {
	last if $line =~ /^$/;
	chomp $line;
	@line = split(/\t/, $line, $num_col);
	$sum_cov = 0;
	$mean_cov = 0;
	$geno_calls = 0;
	$allele_count = 0;
	# foreach individual
	foreach my $depth ( @line[$first_sample_col..$last_sample_col] ){
		if ($depth ne "") { 
			$geno_calls++;
			$allele_count++;
			# if the individual is heterozygous
			if ($depth =~ /\//) { 
				# the following 4 lines will count all reported allele depths
				# even the allele depth of the third and fourth "allele", etc.
				# in a gentoype
				my @allele = split(/\//, $depth);
				foreach my $allele ( @allele ) {
					$sum_cov += $allele;
				}
				$allele_count++;
				# the next three lines will only count allele depth in heterozygous genotype calls
#				$depth =~ /(\d+)\/(\d+)/;	
#				$sum_cov += $1;
#				$sum_cov += $2;
#				die "het allele depth not stored\n" if ( $1 eq "" || $2 eq ""); 
			}else {
				$sum_cov += $depth;
			}
		}
	}
	unless ($allele_count) { 
#		print STDERR "Found no allele depth for catalog locus $line[0]!\n"; next; 
		print ERR join("	", @line[0..$first_sample_col-1], 
			$geno_calls, 
			$mean_cov, 
			"NA",
			@line[$first_sample_col..$#line]
		);
		print ERR "\n";
		next;
	}
	# calculate average read coverage per allele for this locus
#	$mean_cov = $sum_cov/$allele_count;
	# print out
	print join("	", @line[0..$first_sample_col-1], 
			$geno_calls, 
			$mean_cov, 
			$sum_cov,
			@line[$first_sample_col..$#line]
		);
	print "\n";
}
close IN;
close ERR;
