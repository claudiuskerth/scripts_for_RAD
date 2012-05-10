#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: filter_credible_tags.pl
#
#        USAGE: ./filter_credible_tags.pl  <export_output_file_from_database> 
#
#  DESCRIPTION: 
#		This script filters credible SNPs form stacks' export_sql.pl output. 
# 		It filters for a certain value in the deleveraged column and also for only 
# 		up to 2 alleles within a "genotype" assigned by the ML algorithm of stacks.  
# 		More than 2 alleles can be assigned if more than two stacks match one stack 
# 		in the catalogue. Those "genotype" calls are ambiguous and indicate duplicated 
# 		loci.
# 		The script currently expects a sample from two populations and filters for at 
# 		least one genotype call for each population. It also looks for negative Fis values 
# 		within each population as well as in the total sample. Currently, all loci wtih Fis 
# 		below -0.1 will be discarded. This last filter is taken from Hohenlohe2011.
#		The script prints the rows that pass its filters to STDOUT.
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 09/05/12 14:09:07
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

#use Data::Dumper;

my $usage = "
This script is designed to filter ***credible*** loci from the output of \"export_sql.pl\", 
the database export programme of the RAD analysis pipeline \"stacks\". It expects two populations
in the output and all samples from one population must appear left of all samples from the other
population in the export output table. It expects the sample columns to start after the \"Deleveraged\"
column.
The programme filters out loci:

- for which the deleveraging algorithm has been turned on in ustacks for any individual (this is conservative)
- for which the number of catalog matches is higher than the total number of genotype calls. Sometimes two or more stacks
  from an individual can match the same catalog stack. Then the genotype reported for this individual shows alleles from
  both stacks. This is generally an indication of similar duplicated loci having been merged into one catalog stack.
- which have no genotype call in either population
- for which the within population and total sample inbreeding coefficient is < -0.1. Negative F_is values indicate an excess of observed
  heterozygotes over HW expectations, which can be caused by clustering similar loci together.

The script expects a .tsv file and prints input rows which pass the filters to standard output. So pipe it into a file to safe it.

\n$0 <export_output_from_stacks_database>\n
";

die $usage unless @ARGV == 1;
die $usage if ( $ARGV[0] =~ /-h|-help/);

print STDERR "please type in the number of samples from the first population in the export table:\n";
my $pop_a_sample_size = <STDIN>;
chomp $pop_a_sample_size;
die "try again and specify a NUMBER of samples!\n" if $pop_a_sample_size =~ /\D/;

# read in the export file name from the command lines
my $infile = shift;
chomp $infile;
open(IN, "<", "$infile") or die "Can't open the specified file";

# get the first line (header line)
my $line = <IN>;
chomp $line;
my @header = split( /\t/, $line );

# get the index number of the "Alleles" column
my ($alleles_col) = grep { $header[$_] eq "Alleles" } 0..$#header;

# get the index number of the "Deleveraged" column
my ($deleveraged_col) = grep { $header[$_] eq "Deleveraged" } 0..$#header;

# get the index number of the "Num Parents" column
# the "Num Parents" column reports the number of matches from individual
# stacks to the catalog locus
my ($num_pare) = grep { $header[$_] eq "Num Parents" } 0..$#header;

# get the number of samples in the output
my @num_samples =  @header[$deleveraged_col+1..$#header];

# print out headers
print join('	', 
	@header[0..$deleveraged_col],
    "Num_geno_calls",
	"Geno_calls_PopLeft",
	"Geno_calls_PopRight",
	@header[$deleveraged_col+1..$#header],	
	"H_obs_overall", 
	"H_exp_overall", 
	"F_is_overall", 
	"H_obs_PopLeft", 
	"H_exp_PopLeft", 
	"F_is_PopLeft", 
	"H_obs_PopRight", 
	"H_exp_PopRight", 
	"F_is_PopRight"), 
"\n"; 

# initalize global variables
my $delev_count = 0;
my $locus_ok = 0;
my $disploid_count = 0;
my $cat_match = 0; # counts loci for which there are more catalog matches than genotype calls
my $no_PopRight_geno = 0; # counts the number of loci for which there is no genotype call for the right populations in the sql export file
my $no_PopLeft_geno = 0;  # counts the number of loci for which there is no genotype call for the left populations in the sql export file
my $message = "";
my $dodgy_tag_count = 0;
my $consensus_het = 0; # counts loci which contain "consensus/consensus" genotype calls, which is dodgy
 # :TODO      :09/05/12 15:07:08:CEK: insert subroutine, that checks that the number of genotype calls equals the number of parental matches


while (<IN>) {
	# stop reading in the file, when meeting a blank line
	# and write report to standard error 
	if (/^$/) {
		last;
	}else {
		chomp;
		# I have to give the 'split' function a LIMIT here (the length of the header line), 
		# otherwise it would truncate empty fields at the end of the row
		my @row = split( /\t/, $_, @header );
		
	my ($message, $num_geno_calls, $geno_calls_PopLeft, $geno_calls_PopRight) = CHECK_LOCUS(\@row);
		
		if ( $message ) {
			print STDERR $row[0], "\t$message";
		}else {
			$locus_ok++;
			print join("	", @row[0..$deleveraged_col],
			 					$num_geno_calls, 
								$geno_calls_PopLeft, 
								$geno_calls_PopRight,
								@row[$deleveraged_col+1..$#row]
			 		  );
			print "\n";
		}
	}
}

print STDERR "I discarded ", $. - $locus_ok -2, " loci and wrote $locus_ok catalog loci to standard output.\n";
print STDERR "I found ", $delev_count, " deleveraged loci.\n";
print STDERR "I found ", $disploid_count, " loci which had disploid genotype calls.\n";
print STDERR "I found ", $cat_match, " loci which had more catalog matches than genotype calls.\n";
print STDERR "I found ", $no_PopRight_geno, " loci which had no genotype call for PopRight.\n";
print STDERR "I found ", $no_PopLeft_geno, " loci which had no genotype call for PopLeft.\n";
print STDERR "I found ", $dodgy_tag_count, " dodgy tags with F_is < -0.1.\n";
print STDERR "I found ", $consensus_het, " tags with a \"consensus/consensus\" genotype call.\n"; 
close IN;

exit;

##############################################################################
####
#### SUBROUTINES
####
##############################################################################

#===  FUNCTION  ================================================================
#         NAME: CHECK_LOCUS
#      PURPOSE: Checks the catalog locus for filter criteria
#   PARAMETERS: one line from input, as array reference
#      RETURNS: filter message or nothing, if locus passes the filters 
#  DESCRIPTION: 
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
	
sub CHECK_LOCUS {
	
	my $row = shift;

	# check whether there are more matches to the catalog than
	# genotype calls:
	my ($filter_message, $num_geno_calls) = CHECK_CATALOG_MATCHES($row); 
	return $filter_message if $filter_message;

	my %hap_count = (); # stores counts for each haplotype in each pop 

	# store the haplotypes in the whole sample from the "Alleles" column
	my @haplotypes_in_pop = split ';', ${$row}[$alleles_col];

	# if the locus is monomorphic, then the "Alleles" column is empty 
	unless (@haplotypes_in_pop) { push @haplotypes_in_pop, "consensus"; }

	# initialize the hash containing the allele counts
	foreach my $hap (@haplotypes_in_pop) {
		$hap_count{$hap}->{PopLeft} = 0;
		$hap_count{$hap}->{PopRight} = 0; 
	}
	
	# initialize hash for genotype counts
	my %geno_count = (); # stores counts for het and hom genotypes in the sample
	$geno_count{PopRight}->{all} = 0;
	$geno_count{PopRight}->{het} = 0;
	$geno_count{PopRight}->{hom} = 0;
	$geno_count{PopLeft}->{all} = 0;
	$geno_count{PopLeft}->{het} = 0;
	$geno_count{PopLeft}->{hom} = 0;
	
	# take only the loci for which NO deleveraging algorithm 
	# has been turned on in any individual by ustacks 
	if ( ${$row}[$deleveraged_col] == 0  ) {

		$filter_message = FILTER_DISPLOID_GENO(
			{		
				line =>	$row, 
				count_haplotypes => \%hap_count, 
				count_genotypes => \%geno_count,	
			}
		);
#		print Dumper(%geno_count);
		
		# if no filter message has been returned from the previous subroutine
		# and there is at least one gentoype call for each of PopRight and PopLeft, 
		# filter out loci with  negative F_is, i. e. an excess of observed heterozygote
		# genotypes over HW expectations:
		if ( !$filter_message and 
			$geno_count{PopRight}->{all} and $geno_count{PopLeft}->{all} ) {
			$filter_message = FILTER_NEG_FIS(
				{
					count_for_genos => \%geno_count,	
					count_for_haps => \%hap_count,
					line => $row
				}
			);
		# if there is no filter message but also no genotype call for PopRight
		}elsif ( !$filter_message and !$geno_count{PopRight}->{all} ) {
			$filter_message = "No genotype called for PopRight at catalog locus ". ${$row}[0] . ".\n";
			$no_PopRight_geno++;
#			print Dumper(%geno_count), "\n";
		# if there is no filter message but also no genotype call for PopLeft
		}elsif ( !$filter_message and !$geno_count{PopLeft}->{all} ) {
			$filter_message = "No genotype called for PopLeft at catalog locus ". ${$row}[0] . ".\n";
			$no_PopLeft_geno++;
		}
		return($filter_message, $num_geno_calls, $geno_count{PopLeft}->{all}, $geno_count{PopRight}->{all}); 	
		
	# if the deleveraging algorithm has been turned on by ustacks
	}elsif ( ${$row}[$deleveraged_col] > 0 ) { 
		$delev_count++;
		return("Found deleveraged locus ... going to ignore\n");
	}
}

#===  FUNCTION  ================================================================
#         NAME: CHECK_CATALOG_MATCHES
#      PURPOSE: filter out loci where two or more separate individual stacks match
#      			the same catalog locus, which should indicate improper read clustering,
#      			probably due to similar duplicated loci.
#   PARAMETERS: reference to input line
#      RETURNS: nothing or filter message 
#  DESCRIPTION: Count the number of genotype calls in the total sample for this locus and compare it 
# 				to the number of matches to the catalog. If read clustering within individuals and
# 				across individuals for the catalog of loci has worked, then both numbers should be the same.
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================

sub CHECK_CATALOG_MATCHES {

	my	( $row )	= @_;

	my $message = "";
	my $geno_calls = 0;
	
	foreach my $ind ( @{$row}[$deleveraged_col+1..$#{$row}] ) {
		$geno_calls++ if $ind ne "";
	}
	
	if ( $geno_calls < ${$row}[$num_pare] ) {
		$cat_match++;
		$message = "The catalog locus with the ID ${$row}[0] has more matches to the catalog than genotype calls ... going to ignore\n";
		return($message, $geno_calls);
	}
	return($message, $geno_calls);
	
} ## --- end sub CHECK_CATALOG_MATCHES


#===  FUNCTION  ================================================================
#         NAME: FILTER_NEG_FIS
#      PURPOSE: checks the catalog locus for negative Fis values
#   PARAMETERS: genotype count and haplotype count hash as well as input line, all
#   			as references 
#      RETURNS: nothing if locus has 
#  DESCRIPTION: calculates expected and observed heterozygosity within each population
#  				and in the total sample as well as the individual inbreeding coefficient.
#  				The locus is discarded if any of the Fis values is below -0.1, which a
#  				somewhat arbitrary threshold. 
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================

sub FILTER_NEG_FIS {
	
	my $arg_ref = shift;
	
	my $geno_count = $arg_ref->{count_for_genos};
	my $hap_count = $arg_ref->{count_for_haps};
	my $row = $arg_ref->{line};
	
	my $all_genotype_calls = ${$geno_count}{PopLeft}->{all} 
								+ ${$geno_count}{PopRight}->{all};

# calculate OBSERVED heterozygosity for the locus
            
    # overall

    my $H_obs_overall = ( ${$geno_count}{PopRight}->{het} 
     					+ ${$geno_count}{PopLeft}->{het} )
     					/ ( $all_genotype_calls );

	# for PopLeft pop
	my $H_obs_PopLeft= ${$geno_count}{PopLeft}->{het}
							/ ${$geno_count}{PopLeft}->{all};

	# for PopRight pop 
	my $H_obs_PopRight= ${$geno_count}{PopRight}->{het}
					/ ${$geno_count}{PopRight}->{all};
    					
# calculate EXPECTED heterozygosity for the locus
    my $ss_overall = 0;
    my $ss_PopLeft = 0;
    my $ss_PopRight = 0;
    my $Fis_overall = 0;
    my $Fis_PopLeft = 0;
    my $Fis_PopRight = 0;                     
    my $H_exp_overall = 0;
    my $H_exp_PopLeft = 0;
    my $H_exp_PopRight = 0;
    
    # leave exp. Het. and F_is = 0 unless the tag is variable
    unless ( exists ${$hap_count}{consensus} ) {
    	
        foreach my $haplotype ( keys %{$hap_count} ) {
		    # get sum of squares of allele counts
		    $ss_overall += ( ${$hap_count}{$haplotype}->{PopRight} 
		                     + ${$hap_count}{$haplotype}->{PopLeft} )**2;
		
		    $ss_PopLeft += ${$hap_count}{$haplotype}->{PopLeft}**2;
		                                            
		    $ss_PopRight += ${$hap_count}{$haplotype}->{PopRight}**2;
		}
          
    	# see Hedrick, p. 98, example 2.10 for an explanation of the small sample 
    	# size correction for the expected heterozygosity
    	$H_exp_overall = ( $all_genotype_calls*2 ) / ( $all_genotype_calls*2 -1 ) 
    	                   * ( 1 - ( $ss_overall / ($all_genotype_calls*2)**2 ) ); # 1 minus the expected homozygosity

    	$H_exp_PopLeft= ( ${$geno_count}{PopLeft}->{all}*2 ) / ( ${$geno_count}{PopLeft}->{all}*2 -1 ) 
    	        		           * ( 1 - ( $ss_PopLeft / (${$geno_count}{PopLeft}->{all}*2)**2 ) );

		$H_exp_PopRight = ( ${$geno_count}{PopRight}->{all}*2 ) / ( ${$geno_count}{PopRight}->{all}*2 -1 ) 
    	        		           * ( 1 - ( $ss_PopRight / (${$geno_count}{PopRight}->{all}*2)**2 ) );

    }#close unless

	## calculate F_is for the locus

	if ( $H_exp_overall ) { 
	    $Fis_overall = 1-($H_obs_overall/$H_exp_overall);
	    if ( abs( 0 + $Fis_overall) < 0.00000001) { $Fis_overall = 0 ;} 
	}
	if ( $H_exp_PopLeft ) { 
		$Fis_PopLeft = 1-($H_obs_PopLeft/$H_exp_PopLeft); 
		if ( abs( 0 + $Fis_PopLeft ) < 0.00000001) { $Fis_PopLeft = 0 ;}
	}       
	if ( $H_exp_PopRight ) { 
		$Fis_PopRight  = 1-($H_obs_PopRight/$H_exp_PopRight); 
		if ( abs( 0 + $Fis_PopRight ) < 0.00000001) { $Fis_PopRight = 0 ;}
	}                         


	if ( $Fis_overall < -0.1 || $Fis_PopLeft < -0.1 || $Fis_PopRight < -0.1 ) {
	    $dodgy_tag_count++;
	    return("The tag with catalog ID ". ${$row}[0] . " is dodgy!\n");
	}else {
	    push @{$row}, ($H_obs_overall, $H_exp_overall, $Fis_overall, 
	                            $H_obs_PopLeft, $H_exp_PopLeft, $Fis_PopLeft, 
	                            $H_obs_PopRight, $H_exp_PopRight, $Fis_PopRight);
	}   

	return;
} 

#===  FUNCTION  ================================================================
#         NAME: FILTER_DISPLOID_GENO
#      PURPOSE: weeds out loci for which there exist disploid genotype calls 
#   PARAMETERS: input line as array reference, references to haplotype count
#   			and genotype count hashes
#      RETURNS: filter message or nothing if the locus passes the filter
#  DESCRIPTION: iterates over each individual genotype and checks for disploid 
#  				genotype calls, i. e. more than 2 alleles by calling the subroutine
#  				GENOTYPE for each genotype call
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================

sub FILTER_DISPLOID_GENO {
	my $arg_ref = shift;
	
	my $row = $arg_ref->{line};
	my $hap_count = $arg_ref->{count_haplotypes};
	my $geno_count = $arg_ref->{count_genotypes};
	
	# iterate over each individual genotype field
	my $pop = "";
	my $field_count = 0;
	my $message;
		
	foreach my $genotype ( @{$row}[$deleveraged_col+1..$#{$row}] ) {
		$field_count++;
		if ( $genotype ne "" ) {
			
			if ( $field_count <= $pop_a_sample_size ) {
				$pop = "PopLeft";
				${$geno_count}{PopLeft}->{all}++;
			}else { 
				$pop = "PopRight";
				${$geno_count}{PopRight}->{all}++; 
			}
			
			$message = GENOTYPE( 
				{
					genotype => $genotype, 
					population => $pop,
					geno_count => $geno_count,
					hap_count => $hap_count,
					row => $row,
				}
			);
			# If a filter message has been returned for the genotype, stop 
			# iterating over the rest of the genotypes and return that error 
			# message.
			return($message) if ($message);	
		}			
	}# end GENOTYPES
	return($message);	
}

#===  FUNCTION  ================================================================
#         NAME: GENOTYPE
#      PURPOSE: detect genotype calls with more than 2 alleles 
#   PARAMETERS: gets the genotype, population, input line, genotype count and
#   			haplotype count hashes, all as references 
#      RETURNS: filter message or nothing if the genotype is ok 
#  DESCRIPTION: checks the genotype for disploid genotype calls by counting the
#  				number of alleles in the genotype call
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================

sub GENOTYPE {
	my $arg_ref = shift;
	
	my $genotype = $arg_ref->{genotype};
	my $pop = $arg_ref->{population};
	my $geno_count = $arg_ref->{geno_count};
	my $hap_count = $arg_ref->{hap_count};
	my $row = $arg_ref->{row};

	
	# count haplotypes per "genotype" assigned to each indidvidual
	my $num_haplotypes = 1;
	# count the number of "/" in each field
	$num_haplotypes += ($genotype =~ tr[\/][\/]);
	
	# initialize variable that is incremented when a haplotype from 
	# among the list of haplotypes in the sample has been found in 
	# this genotype (this is just to check the proper function of the code)
	my $found = 0;
	
	# filter out catalog loci for which there are triploid or higher 
	# ploid inividual "genotype" calls:
	if ( $num_haplotypes > 2 ) {		
		my $message = "The catalog locus " . 
					  ${$row}[0] .
					  " is associated with disploid individual genotype calls " .
					  "... going to ignore this whole locus\n";
		$disploid_count++;
		return ($message);
	}
	elsif ( $num_haplotypes == 2 ) { 
		if ( exists ${ $hap_count }{consensus} ) { 
			my $message = "The catalog locus " . 
							${$row}[0] .
						  " has a consensus heterozygote genotype call. " .
					  	  "... most  likely dodgy ... going to ignore.\n";	  
			$consensus_het++;
			return ($message);
		}
		# foreach haplotype found in the sample
		else {
			foreach my $hap ( keys %{ $hap_count }) {
				if ( $genotype =~ /$hap/ ) {
					# count the haplotype found once 
					${$hap_count}{$hap}->{$pop}++;
					$found++;
				}
			}
			${$geno_count}{$pop}->{het}++;
		}
	}else { 
		${$geno_count}{$pop}->{hom}++; 
		# foreach haplotype found in the population
		foreach my $hap ( keys %{ $hap_count }) {
			if ( $genotype =~ /$hap/ ) {
				# count the haplotype found twice 
				${$hap_count}{$hap}->{$pop} += 2;
				$found++;
			}
		}
	}
	unless ($found) { 
		print STDERR "Error: Found no haplotype for genotype in subroutine genotype\n" .
				"\tcatalog_no: " . $row->[0] . "\n".
				"\talleles in sample: ". $row->[10] ."\n";
		foreach my $hap ( keys %{ $hap_count }) {		
				print STDERR $hap, "\t";
		}
		print STDERR "\n";
	}
	return;
}# end sub genotype
################################################################################
