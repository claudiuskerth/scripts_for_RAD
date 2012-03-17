#!/usr/bin/env perl
use strict; use warnings;

# get the file name to read
my $infile = shift or die;

# open file to read
open(IN, "<$infile") or die;

# read in the first line
my $line = <IN>;

my @line = split(/\t/, $line);
my $num_col = @line;
my $num_ind = @line[12..$#line];

print join("	", @line[0..11]), "\t", "geno_calls", "\t", "average coverage/allele", "\t";
print join("	", @line[12..$#line]);

my ($sum_cov, $mean_cov, $geno_calls, $allele_count);

while ($line = <IN>) {
	@line = split(/\t/, $line, $num_col);
	$sum_cov = 0;
	$mean_cov = 0;
	$geno_calls = 0;
	$allele_count = 0;
	foreach my $depth (@line[12..$#line]){
		chomp $depth;
		if ($depth ne "") { 
			$geno_calls++;
			$allele_count++;
			if ($depth =~ /\//) { 
				$allele_count++;
				my ($allele1, $allele2) = $depth =~ /(\d+)\/(\d+)/;	
				$sum_cov += $allele1;
				$sum_cov += $allele2;
				die "het allele depth not stores\n" if ( $allele1 eq "" || $allele2 eq ""); 
			}else {
				$sum_cov += $depth;
			}
		} else { next; }
		
	}
	$mean_cov = $sum_cov/$allele_count;
	print join("	", @line[0..11]), "\t", $geno_calls, "\t", $mean_cov, "\t", join("	", @line[12..$#line]), "\n";
}
close IN;












































