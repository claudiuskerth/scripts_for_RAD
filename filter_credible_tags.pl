#!/usr/bin/env perl
use strict; use warnings;

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



my $dodgy_tag_count = 0;
my $disploid_locus = 0;
my $locus_ok = 0;

my $infile = shift;
chomp $infile;

open(IN, "<", "$infile") or die "Can't open the specified file";

$_ = <IN>;
chomp;
print join('	', $_ , "H_obs_overall", "H_exp_overall", "F_is_overall", "H_obs_Greixer", "H_exp_Greixer", "F_is_Greixer", "H_obs_Aunat", "H_exp_Aunat", "F_is_Aunat"), "\n"; 

while (<IN>) {
	if (/^$/) {
		print STDERR "I discarded ", $. - $locus_ok -2, " loci and wrote $locus_ok catalog loci to output.\n";
		last;
	}
	chomp;
	# I have to give the 'split' function a LIMIT here, 
	# otherwise it would truncate empty fields at the end of the row
	my @row = split( /\t/, $_, 48 );
	
	# take only the loci for which NO deleveraging algorithm 
	# has been turned on in any individual
	if ($row[11] == 0 ) {
		# filter out loci for which individuals have tri-, tetra-, etc. ploid 
		# genotype calls, which should indicate ambiguous assignment of alleles 
		# to loci
		my $non_diploid_genotype = 0;
		my $num_hetero_geno = 0;
		my $num_hetero_geno_Greixer = 0;
		my $num_hetero_geno_Aunat = 0;
		my $genotype_calls_Greixer = 0;
		my $genotype_calls_Aunat = 0;
		my $all_genotype_calls = 0;
		my %haplotypes_in_pop;
		my $field_count = 0;
	
#		store the haplotypes in the whole sample from the "Alleles" column
		my @haplotypes_in_pop = split ';', $row[10];

#		initialize the hash
		foreach my $haplotype (@haplotypes_in_pop) {
			$haplotypes_in_pop{$haplotype}->{Greixer} = 0;
			$haplotypes_in_pop{$haplotype}->{Aunat} = 0; 
		}
		
		# iterate over each individual genotype field
		foreach my $field ( @row[12..$#row] ) { 

			$field_count++;
			
			if ( $field ne "" ) { ###						
				$all_genotype_calls++;
				if ( $field_count <= 18 ) {$genotype_calls_Greixer++;}
				else {$genotype_calls_Aunat++;}
								
# 					count haplotypes per "genotype" assigned to each indidvidual
					my $num_haplotypes = 1;

					# count the number of "/" in each field
					$num_haplotypes += $field =~ tr/\//\//; 

					# count heterozygous genotypes
					if ( $num_haplotypes == 2 ) { 

						$num_hetero_geno++;

						if ( $field_count <= 18 ) {$num_hetero_geno_Greixer++;}
			 			else {$num_hetero_geno_Aunat++;}
			 		}

			 		# record triploid and higher ploid genotypes 
			 		if ( $num_haplotypes > 2 ) { $non_diploid_genotype++; }
			 		
			 		if ( $num_haplotypes > 0 && $num_haplotypes < 3) {
				 		
				 		if ( $field_count <= 18 ) {
				 			foreach my $haplotype (@haplotypes_in_pop) {
				 				if ( $field =~ $haplotype) { 
				 					
				 					# if homozygous genotype
				 					if ( $num_haplotypes == 1 ) {
				 					
				 						# count the haplotype found twice 
				 						$haplotypes_in_pop{$haplotype}->{Greixer} += 2;
				 						# break out of the foreach loop, since a homoz. ind. has 
				 						# only one haplotype at that locus. No need to check for 
				 						# the other haplotypes 
				 						last; 
				 					}
				 					else { $haplotypes_in_pop{$haplotype}->{Greixer}++; }
				 				}
				 			}
				 		}
				 		else {
							foreach my $haplotype (@haplotypes_in_pop) {
				 				if ( $field =~ $haplotype) { 
				 					if ( $num_haplotypes == 1 ) { 
				 						$haplotypes_in_pop{$haplotype}->{Aunat} += 2;														last;
			 						}
				 						else { $haplotypes_in_pop{$haplotype}->{Aunat}++; }
				 					}
			 				}
				 		}
					}
				
			}
		}
		
		
		if ( $non_diploid_genotype == 0 && # only if NO non-diploid genotypes have been found
			$genotype_calls_Aunat > 0 && # only if there is at least one genotype call for Aunat
			$row[5] < 37) { # only if number of individuals matching the catalog is less than 37. It can be higher because stacks 'catalog matches' counts matches of individuals stacks against the catalog. In some individuals, two different stacks can be clustered together as two different alleles from one locus. 

# calculate observed heterozygosity for the locus
		
		# overall
			my $H_obs_overall = $num_hetero_geno/$all_genotype_calls;
		# for Greixer pop	
			my $H_obs_Greixer = $num_hetero_geno_Greixer/$genotype_calls_Greixer;
		# for Aunat pop	
			my $H_obs_Aunat = $num_hetero_geno_Aunat/$genotype_calls_Aunat;

# calculate expected heterozygosity for the locus
			my $ss_overall = 0;
			my $ss_Greixer = 0;
			my $ss_Aunat = 0;
			my $F_is_overall = 0;
			my $F_is_Greixer = 0;
			my $F_is_Aunat = 0;			
			my $H_exp_overall = 0;
			my $H_exp_Greixer = 0;
			my $H_exp_Aunat = 0;
			
			# leave F_is = 0 unless the tag is variable
			if ( $row[7] > 0 ) { 
				# foreach haplotype:			
				foreach my $haplotype ( keys %haplotypes_in_pop ) {

					# get sum of squares of allele counts
					$ss_overall += ( $haplotypes_in_pop{$haplotype}->{Greixer} 
									+ $haplotypes_in_pop{$haplotype}->{Aunat} )**2;
				
					$ss_Greixer += $haplotypes_in_pop{$haplotype}->{Greixer}**2;
										
					$ss_Aunat += $haplotypes_in_pop{$haplotype}->{Aunat}**2;
			
##					print "Greixer", "\t", "Aunat", "\n", 
##					"$haplotype: ", $haplotypes_in_pop{$haplotype}->{Greixer}, "\t", 
##					$haplotypes_in_pop{$haplotype}->{Aunat}, "\n\n"; 
					
				}
###				print STDERR $ss_overall, "\n", $ss_Greixer, "\n", $ss_Aunat, "\n";
###				print STDERR $all_genotype_calls, "\n", $genotype_calls_Greixer, "\n", $genotype_calls_Aunat, "\n";

				
# overall:

# I think I took these small sample size corrections from a paper of Hohenlohe and Julian 	
			
				# small sample size correction:
				$H_exp_overall = ( $all_genotype_calls*2 ) / ( $all_genotype_calls*2 -1 ) 
					* ( 1 - ( $ss_overall / ($all_genotype_calls*2)**2 ) ); # 1 minus the expected homozygosity
					
####				print "\n", "exp. Het.: ", "\t", $H_exp_overall, "\n";
####				print "# of alleles: ", "\t", $all_genotype_calls*2, "\n";
####				print "sums of squares of haplotype count: ", "\t", $ss_overall, "\n";
###						

# for Greixer pop:

				# small sample size correction
				$H_exp_Greixer = ( $genotype_calls_Greixer*2 ) / ( $genotype_calls_Greixer*2 -1 ) 
						* ( 1 - ( $ss_Greixer / ($genotype_calls_Greixer*2)**2 ) ); # 1 minus the expected homozygocity
# for Aunat pop:

				# small sample size correction
				$H_exp_Aunat = ( $genotype_calls_Aunat*2 ) / ( $genotype_calls_Aunat*2 -1 )
						* ( 1 - ( $ss_Aunat / ($genotype_calls_Aunat*2)**2 ) ); # 1 minus the expected homozygocity
	
# calculate F_is for the locus

				# sometimes more than one allele has been found in the superparents, 
				# but the individuals have all the same homozygous genotype. Then the 
				# following formula would try to divide by zero.
				if ( $H_exp_overall ) { 
					$F_is_overall = 1-($H_obs_overall/$H_exp_overall);
					if ( abs( 0 + $F_is_overall) < 0.00000001) { $F_is_overall = 0 ;} 
				} 
						 
				if ($H_exp_Greixer > 0)	{ $F_is_Greixer = 1-($H_obs_Greixer/$H_exp_Greixer); 
					if ( abs( 0 + $F_is_Greixer ) < 0.00000001) { $F_is_Greixer = 0 ;}
				}	
					
				if ($H_exp_Aunat > 0)	{ $F_is_Aunat  = 1-($H_obs_Aunat/$H_exp_Aunat); 
					if ( abs( 0 + $F_is_Aunat ) < 0.00000001) { $F_is_Aunat = 0 ;}
				}	
			
			} # close if ($row[7] > 0)
##			
###			print STDERR $H_exp_overall, "\n", $H_exp_Greixer, "\n", $H_exp_Aunat, "\n";
###			
			my $dodgy_tag = "";
			if ( $F_is_overall < -0.1 || 
				$F_is_Greixer < -0.1 || 
				$F_is_Aunat < -0.1 ) {
					print STDERR "The tag with catalog ID ", $row[0], " is dodgy!\n";
					$dodgy_tag = "dodgy!";
					$dodgy_tag_count++;
			}
	
			if ( $dodgy_tag eq "" ) {
				$locus_ok++;
			 	push @row, ($H_obs_overall, $H_exp_overall, $F_is_overall, 
							$H_obs_Greixer, $H_exp_Greixer, $F_is_Greixer, 
							$H_obs_Aunat, $H_exp_Aunat, $F_is_Aunat);
				print join("	", @row), "\n"; # the join function does not understand "\t"
			}
		}
	} else {
		$disploid_locus++;
		print STDERR "The tag with catalog ID ", $row[0], " is disploid!\n";
		}
	
}
print STDERR "I found ", $dodgy_tag_count, " dodgy tags.\n";
print STDERR "I found ", $disploid_locus, " loci with disploid genotype calls.\n";
close IN;

###################################################################################

