#!/usr/bin/env perl
use strict; use warnings;

#use Data::Dumper;

# FILTER CREDIBLE SNPS FROM STACKS EXPORT version 2
#
# This script filters credible SNPs form stacks' export_sql.pl output. 
# It filters for a certain value in the deleveraged column and also for only 
# up to 2 alleles within a "genotype" assigned by the stacks ML algorithm. 
# More than 2 alleles can be assigned if more than two stacks match one stack 
# in the catalogue (I think). Anyway, those "genotype" calls are ambiguous and 
# indicate duplicated loci. It currently filters for at least one genotype
# call for the Aunat population. It also looks for negative Fis values within
# each population as well as in the total sample. Currently, all loci wtih Fis 
# below -0.1 will be discarded. This last filter is taken from Hohenlohe2011.


my $usage = "\n$0 <export_output_from_stacks_database>\n\n";
die $usage unless @ARGV == 1;

# read in the export file from th command lines
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
my ($num_pare) = grep { $header[$_] eq "Num Parents" } 0..$#header;


# print out headers
print join('	', 
	@header , 
	"H_obs_overall", 
	"H_exp_overall", 
	"F_is_overall", 
	"H_obs_Greixer", 
	"H_exp_Greixer", 
	"F_is_Greixer", 
	"H_obs_Aunat", 
	"H_exp_Aunat", 
	"F_is_Aunat"), 
"\n"; 

# initalize global variables
my $delev_count = 0;
my $locus_ok = 0;
my $disploid_count = 0;
my $cat_match = 0; # counts loci for which there are more than 36 matches to the
# catalog
my $no_Aunat_geno = 0;
my $message = "";
my $dodgy_tag_count = 0;


# read in row by row
while (<IN>) {
	# stop reading in the file, when meeting a blank line
	# and write report to standard output
	if (/^$/) {
		last;
	}else {
		chomp;
		# I have to give the 'split' function a LIMIT here, 
		# otherwise it would truncate empty fields at the end of the row
		my @row = split( /\t/, $_, 48 );
		
		$message = check_locus(\@row);
		
		if ( $message ) {
			print STDERR $message;
		}else {
			$locus_ok++;
			print join("	", @row), "\n";
		}
	}
}

print STDERR "I discarded ", $. - $locus_ok -2, " loci and wrote $locus_ok catalog loci to output.\n";
print STDERR "I found ", $delev_count, " deleveraged loci.\n";
print STDERR "I found ", $disploid_count, " loci which had disploid genotype calls.\n";
print STDERR "I found ", $locus_ok, " loci which where ok.\n";
print STDERR "I found ", $cat_match, " loci which had too many catalog matches.\n";
print STDERR "I found ", $no_Aunat_geno, " loci which had no genotype call for Aunat.\n";
print STDERR "I found ", $dodgy_tag_count, " dodgy tags with F_is < 0.1.\n";

close IN;

exit;

##############################################################################
####
#### SUBROUTINES
####
##############################################################################

##############################################################################
#### Name:       CHECK_LOCUS
#### Function:   Checks the catalog locus for filter criteria
#### Parameters: one line from input,
####			 reference to deleveraged counter
#### Returns:    nothing or filter message
##############################################################################
	
sub check_locus {
	
	my $row = shift;
	
	my %hap_count = (); # stores counts for each haplotype in the sample
	# store the haplotypes in the whole sample from the "Alleles" column
	my @haplotypes_in_pop = split ';', ${$row}[$alleles_col];
	# if the locus is monomorphic, then the "Alleles" column is empty 
	unless (@haplotypes_in_pop) { push @haplotypes_in_pop, "consensus"; }
	# initialize the hash containing the allele counts
	foreach my $hap (@haplotypes_in_pop) {
		$hap_count{$hap}->{Greixer} = 0;
		$hap_count{$hap}->{Aunat} = 0; 
	}
	
	my %geno_count = (); # stores counts for het and hom genotypes in the sample
	$geno_count{Aunat}->{all} = 0;
	$geno_count{Aunat}->{het} = 0;
	$geno_count{Aunat}->{hom} = 0;
	$geno_count{Greixer}->{all} = 0;
	$geno_count{Greixer}->{het} = 0;
	$geno_count{Greixer}->{hom} = 0;
	
	# take only the loci for which NO deleveraging algorithm 
	# has been turned on in any individual and for which there only up to 36
	# stacks from inidividuals matching a catalog stack
	if ( ${$row}[$deleveraged_col] == 0  and ${$row}[$num_pare] < 37 ) {

			my $filter_message = filter_disploid_geno(
				{		
					line =>	$row, 
					count_haplotypes => \%hap_count, 
					count_genotypes => \%geno_count,	
				}
			);
#			print Dumper(%geno_count);
			
			# if no filter message has been returned and there is at least one 
			# gentoype call for Aunat, filter for neg_Fis
			if ( !$filter_message and $geno_count{Aunat}->{all} ) {
				$filter_message = filter_neg_Fis(
					{
						count_for_genos => \%geno_count,	
						count_for_haps => \%hap_count,
						line => $row
					}
				);
			# if there is no filter message but also no genotype call for Aunat
			}elsif ( !$filter_message and !$geno_count{Aunat}->{all} ) {
				$no_Aunat_geno++;
				$filter_message = "No genotype called for Aunat at catalog locus ". ${$row}[0] . ".\n";
#				print Dumper(%geno_count), "\n";
			}
			return($filter_message); 	
			
	# if there are more than 36 matches to the catalog	
	}elsif ( ${$row}[$num_pare] >= 37 ) {
		$cat_match++;
		if ( ${$row}[$deleveraged_col] > 0 ) { $delev_count++; }
		return("Found locus with two many matches to the catalog ... going to ignore.\n");

	}elsif ( ${$row}[$deleveraged_col] > 0 ) { 
		$delev_count++;
		return("Found deleveraged locus ... going to ignore\n");
	}
}
##############################################################################
#### Name:       FILTER_NEG_FIS
#### Function:   Checks the catalog locus for negative Fis values
#### Parameters: one line from input,
####			 reference to genotype count and haplotype count hashes
#### Returns:    nothing or filter message
##############################################################################

sub filter_neg_Fis {
	
	my $arg_ref = shift;
	
	my $geno_count = $arg_ref->{count_for_genos};
	my $hap_count = $arg_ref->{count_for_haps};
	my $row = $arg_ref->{line};
	
	my $all_genotype_calls = ${$geno_count}{Greixer}->{all} 
								+ ${$geno_count}{Aunat}->{all};

# calculate OBSERVED heterozygosity for the locus
            
    # overall

    my $H_obs_overall = ( ${$geno_count}{Aunat}->{het} 
     					+ ${$geno_count}{Greixer}->{het} )
     					/ ( $all_genotype_calls );
	# for Greixer pop       
    my $H_obs_Greixer = ${$geno_count}{Greixer}->{het}
    					/ ${$geno_count}{Greixer}->{all};
    # for Aunat pop 
    my $H_obs_Aunat = ${$geno_count}{Aunat}->{het}
    					/ ${$geno_count}{Aunat}->{all};
    					
# calculate EXPECTED heterozygosity for the locus
    my $ss_overall = 0;
    my $ss_Greixer = 0;
    my $ss_Aunat = 0;
    my $Fis_overall = 0;
    my $Fis_Greixer = 0;
    my $Fis_Aunat = 0;                     
    my $H_exp_overall = 0;
    my $H_exp_Greixer = 0;
    my $H_exp_Aunat = 0;
    
    # leave F_is = 0 unless the tag is variable
    unless ( exists ${$hap_count}{consensus} ) {
    	
        foreach my $haplotype ( keys %{$hap_count} ) {
		    # get sum of squares of allele counts
		    $ss_overall += ( ${$hap_count}{$haplotype}->{Greixer} 
		                     + ${$hap_count}{$haplotype}->{Aunat} )**2;
		
		    $ss_Greixer += ${$hap_count}{$haplotype}->{Greixer}**2;
		                                            
		    $ss_Aunat += ${$hap_count}{$haplotype}->{Aunat}**2;
		}
    }#close unless
          
    # see Hedrick, p. 98, example 2.10 for an explanation of the small sample 
    # size correction for the expected heterozygosity
    $H_exp_overall = ( $all_genotype_calls*2 ) / ( $all_genotype_calls*2 -1 ) 
                       * ( 1 - ( $ss_overall / ($all_genotype_calls*2)**2 ) ); 
                       # 1 minus the expected homozygosity
    $H_exp_Greixer = ( ${$geno_count}{Greixer}->{all}*2 ) / ( ${$geno_count}{Greixer}->{all}*2 -1 ) 
                       * ( 1 - ( $ss_Greixer / (${$geno_count}{Greixer}->{all}*2)**2 ) );
                       
    $H_exp_Aunat = ( ${$geno_count}{Aunat}->{all}*2 ) / ( ${$geno_count}{Aunat}->{all}*2 -1 ) 
                       * ( 1 - ( $ss_Aunat / (${$geno_count}{Aunat}->{all}*2)**2 ) );
                       
## calculate F_is for the locus

    # sometimes more than one allele has been found in the superparents, 
    # but the individuals have all the same homozygous genotype. Then the 
    # following formula would try to divide by zero.
    if ( $H_exp_overall ) { 
        $Fis_overall = 1-($H_obs_overall/$H_exp_overall);
        if ( abs( 0 + $Fis_overall) < 0.00000001) { $Fis_overall = 0 ;} 
    }
	if ($H_exp_Greixer > 0) { 
		$Fis_Greixer = 1-($H_obs_Greixer/$H_exp_Greixer); 
		if ( abs( 0 + $Fis_Greixer ) < 0.00000001) { $Fis_Greixer = 0 ;}
    }       
	if ($H_exp_Aunat > 0) { 
		$Fis_Aunat  = 1-($H_obs_Aunat/$H_exp_Aunat); 
		if ( abs( 0 + $Fis_Aunat ) < 0.00000001) { $Fis_Aunat = 0 ;}
	}                         

    if ( $Fis_overall < -0.1 || $Fis_Greixer < -0.1 || $Fis_Aunat < -0.1 ) {
        $dodgy_tag_count++;
        return("The tag with catalog ID ". ${$row}[0] . " is dodgy!\n");
    }else {
        push @{$row}, ($H_obs_overall, $H_exp_overall, $Fis_overall, 
                                $H_obs_Greixer, $H_exp_Greixer, $Fis_Greixer, 
                                $H_obs_Aunat, $H_exp_Aunat, $Fis_Aunat);
   	}   
	return;
} 

##############################################################################
#### Name:       FILTER_DISPLOID_GENO
#### Function:   iterates over each individual genotypes
####			 and checks for disploid genotype calls
#### Parameters: reference to array of row elements,
#### Returns:    nothing or filter message
##############################################################################

sub filter_disploid_geno {
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
			
			if ( $field_count <= 18 ) {
				$pop = "Greixer";
				${$geno_count}{Greixer}->{all}++;
			}else { 
				$pop = "Aunat";
				${$geno_count}{Aunat}->{all}++; 
			}
			
			$message = genotype( 
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

##############################################################################
#### Name:       GENOTYPE
#### Function:   checks genotype for disploid genotype calls
#### Parameters: reference to anonymous hash of arguments,
#### Returns:    nothing or filter message
##############################################################################
sub genotype {
	my $arg_ref = shift;
	
	my $genotype = $arg_ref->{genotype};
	my $pop = $arg_ref->{population};
	my $geno_count = $arg_ref->{geno_count};
	my $hap_count = $arg_ref->{hap_count};
	my $row = $arg_ref->{row};

	
	# count haplotypes per "genotype" assigned to each indidvidual
	my $num_haplotypes = 1;
	# count the number of "/" in each field
	$num_haplotypes += ($genotype =~ tr/\//\//);
	
	# initialize variable that is incremented when a haplotype from 
	# among the list of haplotypes in the sample has been found in 
	# this genotype (this is just to check the proper function of the code)
	my $found = 0;
	
	# filter out catalog loci for which there are triploid or higher 
	# ploid inividual "genotype" calls (ML):
	if ( $num_haplotypes > 2 ) {		
		my $message = "The catalog locus " . 
					  ${$row}[0] .
					  " is associated with disploid individual genotype calls " .
					  "... going to ignore this whole locus\n";
		$disploid_count++;
		return ($message);
	}
	elsif ( $num_haplotypes == 2 ) { 
		${$geno_count}{$pop}->{het}++;
		# foreach haplotype found in the sample
		foreach my $hap ( keys %{ $hap_count }) {
			if ( $genotype =~ /$hap/ ) {
				# count the haplotype found once 
				${$hap_count}{$hap}->{$pop}++;
				$found++;
			}
		}
	}else { 
		${$geno_count}{$pop}->{hom}++; 
		# foreach haplotype found in the population
		foreach my $hap ( keys %{ $hap_count }) {
			if ( $genotype =~ /^$hap$/ ) {
				# count the haplotype found twice 
				${$hap_count}{$hap}->{$pop} += 2;
				$found++;
			}
		}
	}
	unless ($found) { 
		print STDERR "Error: Found no haplotype for genotype in subroutine genotype\n" .
				"\tcatalog_no: " . $row->[0] . "\n";
	}
	return;
}# end sub genotype
################################################################################
