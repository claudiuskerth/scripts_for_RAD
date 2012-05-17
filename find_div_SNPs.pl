#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: find_div_SNPs.pl
#
#        USAGE: ./find_div_SNPs.pl  
#
#  DESCRIPTION: This is script is designed to calculate Weir and Cockerham F_st 
#               values for bi-allelic SNP's in the exported RAD tags from a stacks
#               database. It is certainly going to change in the near future to
#               calculate Nst/Kst.
#               Since SNP's within a tag are usually tighly linked, currently the best 
#               use of this script is probably with the "-p" switch, which makes the script
#               only report the F_st of the most divergent SNP in a tag. With the "-t"
#               switch one can also only report the SNP's with an F_st value above a 
#               certain threshold.
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), claudiuskerth[at]gmail.com
# ORGANIZATION: Sheffield University
#      VERSION: 1.0
#      CREATED: 05/11/2012 16:13:08
#     REVISION: ---
#===============================================================================

# :TODO      :05/14/2012 16:03:10:CEK: insert subroutine that calculates Fst based on 
# whole tag sequences instead of just SNP by SNP, i. e. Nst (nucleotide differentiation)
# :TODO      :05/14/2012 16:02:27:CEK: make script able to handle tri-allelic SNP's

use strict;
use warnings;
use Data::Dumper;

my $usage = "
./filter_export_sql_output [options] -i <infile>

The script expects two populations or samples in the stacks sql export table (this should be
your input file) to calculate F_st for. 
It currently only reports F_st values for SNP\'s that are bi-allelic across the total sample. 
The script prints out a line for each SNP that meets the criteria, i. e. more than one output line
per tag is possible unless you have given the \"-p\" switch. The script ignores tags which contain
any tri- or four-allelic SNP. The output file name is the same as the input but with \'Fst\' added 
before the file extension. 

-i    file name of input file (if present in current directory, otherwise whole path)

options:

-t    F_st threshold, i. e. print out only those SNP's with an F_st above this threshold (default is 0)
-p    print out only the SNP with the highest F_st from each tag, if the Fst is above threshold
-m    minimum number of successful genotype calls required for each population before calculating F_st
\n";

# initialize variables taken from the command line
my $Fst_threshold = 0;
my $min_geno_calls = 0;
my ( $infile, $only_most_div_SNP );

PARSE_COMMAND_LINE();

unless ( -e $infile ) { print "the input file you specified doesn't exist\n"; }

open(IN, "< $infile") or die "$!";

# extract input file name and file extension 
my ($infile_stub, $file_ext) = $infile =~ /^(.+)(\..*)$/;
#print $infile_stub, "\n", $file_ext, "\n"; exit;

# open output file
open(OUT, ">", $infile_stub . "_Fst" . $file_ext ) or die "$!";

# get header line from input file
my $header = <IN>;
chomp $header;
my @header = split(/\t/, $header);

## get the index for the column of the leftmost sample in the table 
print STDERR "please type in the name of the leftmost sample in the table (at least its unique part):\n";
my $first_sample_name = <STDIN>;
chomp $first_sample_name;
#my $first_sample_name = "ery_30-10";
my ($first_sample_col) = grep { $header[$_] =~ /$first_sample_name/  } 0..$#header;

## get the index for the column of the rightmost sample in the table
print STDERR "please type in the name of the rightmost sample in the table (at least its unique part):\n";
my $last_sample_name = <STDIN>;
chomp $last_sample_name;
#my $last_sample_name = "par_34-9"; 
my ($last_sample_col) = grep { $header[$_] =~ /$last_sample_name/  } 0..$#header;

# print out new header
print OUT join('	', @header[0..$first_sample_col-1], 
				"SNP_pos", "mean_within_pop_het", "F_st", 
				@header[$first_sample_col..$#header]
		); 
print OUT "\n"; 

## get the sample size for the left population in the table
print STDERR "please type in the number of samples from the first/left population in the export table:\n";
my $PopLeft_sample_size = <STDIN>;
chomp $PopLeft_sample_size;
die "try again and specify a NUMBER of samples!\n" if $PopLeft_sample_size =~ /\D/;
#my $PopLeft_sample_size = 19; 


# get the index for the column that reports the number of SNPs found at the locus
my ($num_snps_col) = grep { $header[$_] =~ /Num SNPs/ } 0..$#header;

my $total_SNP_count = 0;

# read in line by line
LINE:
while(<IN>) {
###	last if $. > 200;
	chomp;
	my @row = split(/\t/, $_);

	# ignore invariable tags
#	next if $row[$num_snps_col] == 0;

	my $num_SNPs = $row[$num_snps_col]; 

	# get the genotypes for each population by taking slices from the input line
	my @PopLeft_geno_calls = @row[$first_sample_col..$first_sample_col+$PopLeft_sample_size-1];
	my @PopRight_geno_calls = @row[$first_sample_col+$PopLeft_sample_size..$last_sample_col];

	# initialize the hash that stores the SNP allele frequencies for each population
	my %SNP_allele_freq	= ();

	my ($geno_calls_PopLeft, $consensus) = SNP_ALLELE_FREQ(
		{
			geno_calls => \@PopLeft_geno_calls,
			num_snps => $num_SNPs,
			snp_allele_freq => \%SNP_allele_freq,
			pop => "PopLeft"
		}
	);

	# if the locus is not variable, do a shortcut
	if ( defined($consensus) ) {
		print OUT join('	', @row[0..$first_sample_col-1],
						"NA", "0", "0", 
						@row[$first_sample_col..$#row]);
		print OUT "\n";				
		next LINE;
	}	

	my ( $geno_calls_PopRight ) = SNP_ALLELE_FREQ(
		{
			geno_calls => \@PopRight_geno_calls,
			num_snps => $num_SNPs,
			snp_allele_freq => \%SNP_allele_freq,
			pop => "PopRight"
		}
	);

	# if there is a minimum of ? genotype calls in either population
	# calculate F_st SNP-wise and print out a separate line for each SNP 
	if ($geno_calls_PopLeft >= $min_geno_calls && $geno_calls_PopRight >= $min_geno_calls){
		CALCULATE_FST(
			{
				num_snps => $num_SNPs,
				snp_allele_freq => \%SNP_allele_freq,
				row => \@row
			}
		);
	}
	else { print "The locus with catalog id ", $row[0], " doesn't have enough genotype calls.\n"; }
}	
close OUT;
exit;
	
	
#===  FUNCTION  ================================================================
#         NAME: SNP_ALLELE_FREQ
#      PURPOSE: calculate allele frequencies for each allele of each SNP in the
#               population
#   PARAMETERS: genotype calls (including empty ones) for the population, the number
#               of SNPs in the tag as reported by stacks, a reference to the 
#               SNP_allele_freq_hash and the working name of the population (for 
#               internal use only)
#      RETURNS: count of successful genotype calls in the population 
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub SNP_ALLELE_FREQ {
	
	my ( $arg_ref ) = @_;

	my $geno_calls = $arg_ref->{geno_calls};
	my $num_SNPs = $arg_ref->{num_snps};
	my $SNP_allele_freq_hash = $arg_ref->{snp_allele_freq};
	my $pop = $arg_ref->{pop};

	# initialise variables
	my @alleles_in_Pop;
	my $geno_calls_count;

	# iterate over each individual genotype from the population
	# and store its alleles in an array
	foreach my $geno ( @{ $geno_calls } ) {
		# export_sql.pl leaves uncalled genotype fields empty, 
		# whereas 'genotypes' gives them a dash '-'. Be sure to 
		# change the next statement when you give the programme 
		# a ML or corrected genotypes input
		unless ($geno eq "") {
			# count successful genotype calls
			$geno_calls_count++; 
			
			# reformat homozygote genotypes from say AT to AT/AT:
			if ( $geno =~ /^([GATC]+)$/ ) { $geno = "$1/$1"; } 			
			elsif ( $geno =~ /consensus/ ) { return($geno_calls_count, "consensus"); }

			# split genotype into it's two alleles
			my @alleles = split(/\//, $geno); 

			# split each allele into an anonymous array and push 
			# it to the array that stores all alleles for this population
			push @alleles_in_Pop, [ split(//, $alleles[0]) ];
			push @alleles_in_Pop, [ split(//, $alleles[1]) ]
		}
	}

	# count SNP alleles for each SNP position:
	#	
	# foreach SNP position
	foreach my $SNP_position ( 0 .. $num_SNPs-1 ) {
		# foreach allele in the population 
		foreach my $allele ( 0 .. $#alleles_in_Pop ) { 
			# count each SNP allele at each SNP position in the population
			${ $SNP_allele_freq_hash }{$pop}->{$SNP_position}->{ $alleles_in_Pop[$allele]->[$SNP_position] }++; 
		}
	}

	# transform SNP allele counts into SNP allele frequencies:
	#
	# foreach SNP position in the tag
	foreach my $SNP_position ( sort keys %{ $SNP_allele_freq_hash->{$pop} } ) {
		# foreach SNP allele at that SNP positiion in the population
		foreach my $SNP_allele ( keys %{ $SNP_allele_freq_hash->{$pop}->{$SNP_position} } ) {
			$SNP_allele_freq_hash->{$pop}->{$SNP_position}->{$SNP_allele} /= ($geno_calls_count * 2);
		}
	} 

	return ($geno_calls_count);
} ## --- end sub SNP_ALLELE_FREQ
		
#===  FUNCTION  ================================================================
#         NAME: CALCULATE_FST
#      PURPOSE: calculates Fst and mean of within population expected heterozygosity
#               SNP-wise for each tag and after checking whether a SNP is bi-allelic.
#   PARAMETERS: the number of SNPs as reported by stacks, a reference to the
#               SNP_allele_freq_hash and a reference to the input line
#      RETURNS: nothing, prints out to file 
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub CALCULATE_FST {
	my	( $arg_ref )	= @_;

	my $num_SNPs = $arg_ref->{num_snps};
	my $SNP_allele_freq = $arg_ref->{snp_allele_freq};
#	print Dumper($SNP_allele_freq);
	my $row = $arg_ref->{row};

	# initialize varibles
	my ($p_Pop_a, $p_Pop_b, @SNP_alleles_in_sample, %SNP_alleles_in_sample, $Fst, $mean_within_pop_het, %Fst_hash); 

	#-------------------------------------------------------------------------------
	#  calculate Fst for each SNP in the tag
	#-------------------------------------------------------------------------------
	# foreach SNP in the tag
	SNP:
   	foreach my $SNP_pos ( 0 .. $num_SNPs-1 ) {
		$total_SNP_count++;	
		# re-initialize, i. e. empty
		%SNP_alleles_in_sample = ();
		@SNP_alleles_in_sample = ();
		# get the SNP_alleles from each population
		foreach my $pop ( keys %{ $SNP_allele_freq } ) {
			foreach my $SNP_allele ( sort keys %{ $SNP_allele_freq->{$pop}->{$SNP_pos} } ) {
				$SNP_alleles_in_sample{$SNP_allele} = 1;
			}
		}
		@SNP_alleles_in_sample = sort keys %SNP_alleles_in_sample;

		# check for tri-allelic SNPs, currently my formula for Fst works only for bi-allelic SNPs
		# i. e. only up to 2 SNP alleles in the total sample 
		if ( @SNP_alleles_in_sample > 2 ) { 
			print "The SNP number ", $SNP_pos+1, " in catalog tag ", ${$row}[0], " has more than two SNP alleles.\n";
			return; # stop iterating over the SNPs of this tag and don't print out Fst for any SNPs of this tag
		} 
		# pick one of the two SNP alleles and if the SNP allele does occur in the population, 
		# then store its allele frequency in a separate variable
		if ( exists $SNP_allele_freq->{PopLeft}->{$SNP_pos}->{ $SNP_alleles_in_sample[0] } ) {
			$p_Pop_a = $SNP_allele_freq->{PopLeft}->{ $SNP_pos}->{ $SNP_alleles_in_sample[0] };
		}else{ $p_Pop_a = 0;}

		if ( exists $SNP_allele_freq->{PopRight}->{$SNP_pos}->{ $SNP_alleles_in_sample[0] } ) {
			$p_Pop_b = $SNP_allele_freq->{PopRight}->{$SNP_pos}->{ $SNP_alleles_in_sample[0] };
		}else{ $p_Pop_b = 0;}

		# calculate Fst and mean with population expected heterozygosity
	    ($Fst, $mean_within_pop_het) = FST(
				{
					p_pop_a => $p_Pop_a,
					p_pop_b => $p_Pop_b,
				}
		);
		# store Fst and mean expected het for each SNP in the tag
		$Fst_hash{$SNP_pos} = [ $Fst, $mean_within_pop_het ]; 

	}# close foreach SNP_pos

	#-------------------------------------------------------------------------------
	#  print out
	#-------------------------------------------------------------------------------
	if ( $only_most_div_SNP ) {
		# sort Fst_hash according to Fst values in descending order and the according to SNP position
		# in ascending order
		foreach my $SNP_pos ( sort { $Fst_hash{$b}[0] <=> $Fst_hash{$a}[0] || $a <=> $b } keys %Fst_hash ) {
			# if the most divergent SNP has an Fst above threshold
			if ( $Fst_hash{$SNP_pos}[0] > $Fst_threshold ) {
				# print out the data only for the SNP with the highest Fst
				# if there 2 or more SNPs with the highest Fst, the most upstream
				# of those is printed out
				print OUT join('	', @{$row}[0..$first_sample_col-1],
								$SNP_pos+1, $Fst_hash{$SNP_pos}[1], $Fst_hash{$SNP_pos}[0],
								@{$row}[$first_sample_col..$#{$row}]);
				print OUT "\n";					
			}
			last; # stop after printing out the SNP with the highest Fst
		}
	}
	else {
		# sort Fst_hash according to Fst values in descending order and the according to SNP position
		# in ascending order
		foreach my $SNP_pos ( sort { $Fst_hash{$b}[0] <=> $Fst_hash{$a}[0] || $a <=> $b } keys %Fst_hash ) {
			# print out all SNPs with an Fst above threshold
			if ( $Fst_hash{$SNP_pos}[0] > $Fst_threshold ) {
				print OUT join('	', @{$row}[0..$first_sample_col-1],
								$SNP_pos+1, $Fst_hash{$SNP_pos}[1], $Fst_hash{$SNP_pos}[0],
								@{$row}[$first_sample_col..$#{$row}]);
				print OUT "\n";					
			}
		}
	}
	return ;
} ## --- end sub CALCULATE_FST


#===  FUNCTION  ================================================================
#         NAME: FST
#      PURPOSE: calculate Fst and mean within population heterozygosity for a SNP
#   PARAMETERS: population SNP allele frequencies
#      RETURNS: Fst and mean within pop het
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub FST {

	my	( $arg_ref )	= @_;

	my $p_Pop_a = $arg_ref->{p_pop_a};
	my $p_Pop_b = $arg_ref->{p_pop_b};

	# initialize $F_st
	my $F_st = 0;
	
	# get the mean of both population SNP allele frequencies
	my $p_mean = ($p_Pop_a+ $p_Pop_b)/2;

	# if the SNP is not variable at all,
	# the Fst is zero, no need to calculate it
	unless ( $p_mean == 1 ) {
		# calculate Weir and Cockerham Fst 
		# I took this formula from the bottom of page 489 in Hedrick's "Genetics of Populations".
		# see also Weir1996, p. 166, for another explanation of this formula 
		$F_st = ( ($p_Pop_a**2 + $p_Pop_b**2)/2 - $p_mean**2 )
				/ ( $p_mean * (1-$p_mean) );

		# the following formula is taken from Weir1996 but it's Fst values are twice as high as the
		# ones calculated by the above formula
#		$F_st = ( ($p_Pop_a - $p_mean)**2 + ($p_Pop_b - $p_mean)**2 )
#				/ ( $p_mean * (1-$p_mean) );
	}
	
#-------------------------------------------------------------------------------
# calculate average within population expected heterozygosity for the SNP
#-------------------------------------------------------------------------------
	# Version 1:
#	my $within_pop_a_het = 2 * $p_Pop_a * (1-$p_Pop_a);	
#	my $within_pop_b_het = 2 * $p_Pop_b * (1-$p_Pop_b);
#	my $mean_within_pop_het = ( $within_pop_a_het + $within_pop_b_het ) / 2;

	# Version 2:
	my $mean_within_pop_het  = 2 * $p_mean * (1-$p_mean);
	
		
	return ( $F_st, $mean_within_pop_het) ;
} ## --- end sub FST


#===  FUNCTION  ================================================================
#         NAME: PARSE_COMMAND_LINE
#      PURPOSE: 
#   PARAMETERS: ????
#      RETURNS: ????
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub PARSE_COMMAND_LINE {
	my $i = 0;
	while(@ARGV){
		$_ = shift @ARGV;
		if ( $_ =~ /^-h|help$/ ) { die $usage; }
		if ( $_ =~ /^-p$/ ) { $only_most_div_SNP = "True"; }
		if ( $_ =~ /^-t$/ ) { $Fst_threshold = shift @ARGV; }
		if ( $_ =~ /^-i$/ ) { $infile = shift @ARGV; $i++; }
		if ( $_ =~ /^-m$/ ) { $min_geno_calls = shift @ARGV; }
	}
	die "You have to specify the input file with the \"-i\" switch.\n$usage" unless $i;
}
