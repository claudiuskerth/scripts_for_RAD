#!/usr/bin/env perl
use strict; use warnings;

my $usage = "
$0 <input_file_with_allele_depths>
Prints to STDOUT, i. e. pipe into file!
";

die $usage unless @ARGV == 1;

# get the file name to read
my $infile = shift or die;

# open file to read
open(IN, "<$infile") or die;

# read in the first line
my $line = <IN>;
my @line = split(/\t/, $line);
my $num_col = @line;

# get the index for the deleveraged column
my ($delev_col) = grep { $line[$_] eq "Deleveraged" } 0..$#line;

my $num_ind = @line[$delev_col+1..$#line];

# print headers and add tow new columns
print join("	", @line[0..$delev_col], 
		"geno_calls", 
		"average coverage/allele", 
		 @line[$delev_col+1..$#line]);

my ($sum_cov, $mean_cov, $geno_calls, $allele_count);
# while reading in loci
while ($line = <IN>) {
	chomp $line;
	@line = split(/\t/, $line, $num_col);
	$sum_cov = 0;
	$mean_cov = 0;
	$geno_calls = 0;
	$allele_count = 0;
	# foreach individual
	foreach my $depth (@line[$delev_col+1..$#line]){
		if ($depth ne "") { 
			$geno_calls++;
			$allele_count++;
			# if the individual is heterozygous
			if ($depth =~ /\//) { 
				$allele_count++;
				$depth =~ /(\d+)\/(\d+)/;	
				$sum_cov += $1;
				$sum_cov += $2;
				die "het allele depth not stored\n" if ( $1 eq "" || $2 eq ""); 
			}else {
				$sum_cov += $depth;
			}
		} else { next; }
		
	}
	# calculate average read coverage per allele for this locus
	$mean_cov = $sum_cov/$allele_count;
	# print out
	print join("	", @line[0..$delev_col], 
			$geno_calls, 
			$mean_cov, 
			@line[$delev_col+1..$#line],
			"\n");
}
close IN;












































