#!/usr/bin/env perl
use strict; use warnings;

# FILTER CREDIBLE SNPS FROM STACKS EXPORT version 2
#
# This version is for parsing corrected genotypes. Empty fields in 'genotypes' 
# output contain a dash, whereas in export_sql.pl output they are actually empty.
# This script filters credible SNPs form stacks' genotypes export output. 
# It filters for a certain value in the deleveraged column and also for only 
# up to 2 alleles within a "genotype" assigned by the stacks ML algorithm. 
# More than 2 alleles can be assigned if more than two stacks match one stack 
# in the catalogue (I think). Anyway, those "genotype" calls are ambiguous and 
# indicate duplicated loci.
#

#my $dodgy_tag_count = 0;
#my $disploid_locus = 0;
#my $locus_ok = 0;
#
#my $infile = shift;
#chomp $infile;
#
#open(IN, "<", "$infile") or die "Can't open the specified file";
#
#$_ = <IN>;
#chomp;
#print join('	', $_ , "H_obs_overall", "H_exp_overall", "F_is_overall", "H_obs_Greixer", "H_exp_Greixer", "F_is_Greixer", "H_obs_Aunat", "H_exp_Aunat", "F_is_Aunat"), "\n"; 
#
#while (<IN>) {
#	if (/^$/) {
#		print STDERR "I discarded ", $. - $locus_ok -2, " loci and wrote $locus_ok catalog loci to output.\n";
#		last;
#	}
#	chomp;
#	# I have to give the 'split' function a LIMIT here, 
#	# otherwise it would truncate empty fields at the end of the row
#	my @row = split( /\t/, $_, 48 );
#	
#	# take only the loci for which NO deleveraging algorithm 
#	# has been turned on in any individual
#	if ($row[11] == 0 ) {
#		# filter out loci for which individuals have tri-, tetra-, etc. ploid 
#		# genotype calls, which should indicate ambiguous assignment of alleles 
#		# to loci
#		my $non_diploid_genotype = 0;
#		my $num_hetero_geno = 0;
#		my $num_hetero_geno_Greixer = 0;
#		my $num_hetero_geno_Aunat = 0;
#		my $genotype_calls_Greixer = 0;
#		my $genotype_calls_Aunat = 0;
#		my $all_genotype_calls = 0;
#		my %haplotypes_in_pop;
#		my $field_count = 0;
#	
##		store the haplotypes in the whole sample from the "Alleles" column
#		my @haplotypes_in_pop = split ';', $row[10];
#
##		initialize the hash
#		foreach my $haplotype (@haplotypes_in_pop) {
#			$haplotypes_in_pop{$haplotype}->{Greixer} = 0;
#			$haplotypes_in_pop{$haplotype}->{Aunat} = 0; 
#		}
#		
#		# iterate over each individual genotype field
#		foreach my $field ( @row[12..$#row] ) { 
#
#			$field_count++;
#			
#			if ( $field ne "" ) { ###						
#				$all_genotype_calls++;
#				if ( $field_count <= 18 ) {$genotype_calls_Greixer++;}
#				else {$genotype_calls_Aunat++;}
#								
## 					count haplotypes per "genotype" assigned to each indidvidual
#					my $num_haplotypes = 1;
#
#					# count the number of "/" in each field
#					$num_haplotypes += $field =~ tr/\//\//; 
#
#					# count heterozygous genotypes
#					if ( $num_haplotypes == 2 ) { 
#
#						$num_hetero_geno++;
#
#						if ( $field_count <= 18 ) {$num_hetero_geno_Greixer++;}
#			 			else {$num_hetero_geno_Aunat++;}
#			 		}
#
#			 		# record triploid and higher ploid genotypes 
#			 		if ( $num_haplotypes > 2 ) { $non_diploid_genotype++; }
#			 		
#			 		if ( $num_haplotypes > 0 && $num_haplotypes < 3) {
#				 		
#				 		if ( $field_count <= 18 ) {
#				 			foreach my $haplotype (@haplotypes_in_pop) {
#				 				if ( $field =~ $haplotype) { 
#				 					
#				 					# if homozygous genotype
#				 					if ( $num_haplotypes == 1 ) {
#				 					
#				 						# count the haplotype found twice 
#				 						$haplotypes_in_pop{$haplotype}->{Greixer} += 2;
#				 						# break out of the foreach loop, since a homoz. ind. has 
#				 						# only one haplotype at that locus. No need to check for 
#				 						# the other haplotypes 
#				 						last; 
#				 					}
#				 					else { $haplotypes_in_pop{$haplotype}->{Greixer}++; }
#				 				}
#				 			}
#				 		}
#				 		else {
#							foreach my $haplotype (@haplotypes_in_pop) {
#				 				if ( $field =~ $haplotype) { 
#				 					if ( $num_haplotypes == 1 ) { 
#				 						$haplotypes_in_pop{$haplotype}->{Aunat} += 2;														last;
#			 						}
#				 						else { $haplotypes_in_pop{$haplotype}->{Aunat}++; }
#				 					}
#			 				}
#				 		}
#					}
#				
#			}
#		}
#		
#		
#		if ( $non_diploid_genotype == 0 && # only if NO non-diploid genotypes have been found
#			$genotype_calls_Aunat > 0 && # only if there is at least one genotype call for Aunat
#			$row[5] < 37) { # only if number of individuals matching the catalog is less than 37. It can be higher because stacks 'catalog matches' counts matches of individuals stacks against the catalog. In some individuals, two different stacks can be clustered together as two different alleles from one locus. 
#
## calculate observed heterozygosity for the locus
#		
#		# overall
#			my $H_obs_overall = $num_hetero_geno/$all_genotype_calls;
#		# for Greixer pop	
#			my $H_obs_Greixer = $num_hetero_geno_Greixer/$genotype_calls_Greixer;
#		# for Aunat pop	
#			my $H_obs_Aunat = $num_hetero_geno_Aunat/$genotype_calls_Aunat;
#
## calculate expected heterozygosity for the locus
#			my $ss_overall = 0;
#			my $ss_Greixer = 0;
#			my $ss_Aunat = 0;
#			my $F_is_overall = 0;
#			my $F_is_Greixer = 0;
#			my $F_is_Aunat = 0;			
#			my $H_exp_overall = 0;
#			my $H_exp_Greixer = 0;
#			my $H_exp_Aunat = 0;
#			
#			# leave F_is = 0 unless the tag is variable
#			if ( $row[7] > 0 ) { 
#				# foreach haplotype:			
#				foreach my $haplotype ( keys %haplotypes_in_pop ) {
#
#					# get sum of squares of allele counts
#					$ss_overall += ( $haplotypes_in_pop{$haplotype}->{Greixer} 
#									+ $haplotypes_in_pop{$haplotype}->{Aunat} )**2;
#				
#					$ss_Greixer += $haplotypes_in_pop{$haplotype}->{Greixer}**2;
#										
#					$ss_Aunat += $haplotypes_in_pop{$haplotype}->{Aunat}**2;
#			
###					print "Greixer", "\t", "Aunat", "\n", 
###					"$haplotype: ", $haplotypes_in_pop{$haplotype}->{Greixer}, "\t", 
###					$haplotypes_in_pop{$haplotype}->{Aunat}, "\n\n"; 
#					
#				}
####				print STDERR $ss_overall, "\n", $ss_Greixer, "\n", $ss_Aunat, "\n";
####				print STDERR $all_genotype_calls, "\n", $genotype_calls_Greixer, "\n", $genotype_calls_Aunat, "\n";
#
#				
## overall:
#
## I think I took these small sample size corrections from a paper of Hohenlohe and Julian 	
#			
#				# small sample size correction:
#				$H_exp_overall = ( $all_genotype_calls*2 ) / ( $all_genotype_calls*2 -1 ) 
#					* ( 1 - ( $ss_overall / ($all_genotype_calls*2)**2 ) ); # 1 minus the expected homozygosity
#					
#####				print "\n", "exp. Het.: ", "\t", $H_exp_overall, "\n";
#####				print "# of alleles: ", "\t", $all_genotype_calls*2, "\n";
#####				print "sums of squares of haplotype count: ", "\t", $ss_overall, "\n";
####						
#
## for Greixer pop:
#
#				# small sample size correction
#				$H_exp_Greixer = ( $genotype_calls_Greixer*2 ) / ( $genotype_calls_Greixer*2 -1 ) 
#						* ( 1 - ( $ss_Greixer / ($genotype_calls_Greixer*2)**2 ) ); # 1 minus the expected homozygocity
## for Aunat pop:
#
#				# small sample size correction
#				$H_exp_Aunat = ( $genotype_calls_Aunat*2 ) / ( $genotype_calls_Aunat*2 -1 )
#						* ( 1 - ( $ss_Aunat / ($genotype_calls_Aunat*2)**2 ) ); # 1 minus the expected homozygocity
#	
## calculate F_is for the locus
#
#				# sometimes more than one allele has been found in the superparents, 
#				# but the individuals have all the same homozygous genotype. Then the 
#				# following formula would try to divide by zero.
#				if ( $H_exp_overall ) { 
#					$F_is_overall = 1-($H_obs_overall/$H_exp_overall);
#					if ( abs( 0 + $F_is_overall) < 0.00000001) { $F_is_overall = 0 ;} 
#				} 
#						 
#				if ($H_exp_Greixer > 0)	{ $F_is_Greixer = 1-($H_obs_Greixer/$H_exp_Greixer); 
#					if ( abs( 0 + $F_is_Greixer ) < 0.00000001) { $F_is_Greixer = 0 ;}
#				}	
#					
#				if ($H_exp_Aunat > 0)	{ $F_is_Aunat  = 1-($H_obs_Aunat/$H_exp_Aunat); 
#					if ( abs( 0 + $F_is_Aunat ) < 0.00000001) { $F_is_Aunat = 0 ;}
#				}	
#			
#			} # close if ($row[7] > 0)
###			
####			print STDERR $H_exp_overall, "\n", $H_exp_Greixer, "\n", $H_exp_Aunat, "\n";
####			
#			my $dodgy_tag = "";
#			if ( $F_is_overall < -0.1 || 
#				$F_is_Greixer < -0.1 || 
#				$F_is_Aunat < -0.1 ) {
#					print STDERR "The tag with catalog ID ", $row[0], " is dodgy!\n";
#					$dodgy_tag = "dodgy!";
#					$dodgy_tag_count++;
#			}
#	
#			if ( $dodgy_tag eq "" ) {
#				$locus_ok++;
#			 	push @row, ($H_obs_overall, $H_exp_overall, $F_is_overall, 
#							$H_obs_Greixer, $H_exp_Greixer, $F_is_Greixer, 
#							$H_obs_Aunat, $H_exp_Aunat, $F_is_Aunat);
#				print join("	", @row), "\n"; # the join function does not understand "\t"
#			}
#		}
#	} else {
#		$disploid_locus++;
#		print STDERR "The tag with catalog ID ", $row[0], " is disploid!\n";
#		}
#	
#}
#print STDERR "I found ", $dodgy_tag_count, " dodgy tags.\n";
#print STDERR "I found ", $disploid_locus, " loci with disploid genotype calls.\n";
#close IN;

###################################################################################

#####################################################################
use Getopt::Std;
our ($opt_t, $opt_h);
getopts('ht:');

# FIND_DIV_SNPS version 2

# This script takes the output from FILTER CREDIBLE SNPS FROM STACKS EXPORT version 2 and filters for
# highly divergent RAD tags between Greixer and Aunat based on allele frequencies of individual SNPs within tags,
# i. e. at least one SNP of one or more SNPs within a tag has to be divergent. I could probably find more divergent
# tags by filtering for divergent tag sequences, i. e. the combination of SNPs within a tag. However, recombination within 
# tags can make distinction between migrants and recombinant hybrids within the hybrid zone uncertain (Ludovic). 


## Still to do:
## repeated code for Aunat and Greixer should be put into subroutines

my $usage = "
usage:
	./filter_export_sql_output [options] infile
		options:
				-t provide the Fst threshold
";

if ($opt_h) { print $usage; exit;}

# get the divergence threshold from the command line
my $threshold = $opt_t;

my $infile = shift;

unless ( -e $infile ) { print "the input file you specified doesn't exist\n"; }

open(IN, "< $infile") or die "$!";

# delete the last 4 characters in the string
my $outfile = substr($infile, 0, -4); # delete the last 4 characters in the string

open(OUT, ">", "$outfile" . "_div_SNPs.tsv" ) or die "$!";


$_ = <IN>;
chomp;
my @columns = split(/\t/, $_);
print OUT join('	', @columns[0..11]), "\t", "genotype calls Aunat", "\t", "genotype calls Greixer", "\t", "div_SNP_pos", "\t", "F_st", "\t",
		join('	', @columns[12..$#columns]), "\n"; 

while(<IN>) {
###	last if $. > 200;
	chomp;
	my @columns = split(/\t/, $_);
		next if $columns[7] == 0;
		my $num_SNPs = $columns[7]; 
###		print scalar @columns, "\n";
#	
		my @Greixer_genotypes = @columns[12..29];
		my @Aunat_genotypes = @columns[30..47];
###		print join("	", @Greixer_genotypes), "\n";
#
		my $genotype_calls_Aunat = 0;
		my $genotype_calls_Greixer = 0;
		my @alleles_Aunat = ();
		my @alleles_Greixer = ();
		my %SNP_allele_freq_Aunat = ();
		my %SNP_allele_freq_Greixer = ();
#		
### Aunat	
		foreach my $ind_genotype (@Aunat_genotypes) {
# export_sql.pl leaves uncalled genotype fields empty, 
# whereas 'genotypes' gives them a dash '-'. Be sure to 
# change the next statement when you give the programme 
# a ML or corrected genotypes input
			if ($ind_genotype ne "") {
				
				$genotype_calls_Aunat++; 
				
				# reformat homozygote genotypes from say AT to AT/AT:
				if ($ind_genotype =~ /^([GATC]+)$/) {$ind_genotype = "$1/$1";} 			
##				print $ind_genotype . "\n";

				# split genotype into it's two alleles
				my @alleles = split(/\//, $ind_genotype); 
##				print "@alleles",  "\n"; 

				foreach my $allele (@alleles) {
					# each allele is split into an anonymous array
					push @alleles_Aunat, [ split(//, $allele) ]; 				}
			}
		}
###			print $haplotypes_Aunat[0]->[0], "\n";
###			print @{$haplotypes_Aunat[0]}, "\n";
###			for my $hap (@haplotypes_Aunat) {
###				print "@{$hap}", "\n";
###			}
#

		# count SNP alleles for each SNP position	
		foreach my $SNP_position ( 0 .. $num_SNPs-1 ) {
			# iterate over the alleles in the Aunat population sample
			foreach my $allele ( 0 .. $#alleles_Aunat ) { 

###				print $allele, " ", $SNP_position, "\n";
###				print $alleles_Aunat[$allele]->[$SNP_position], "\n";

				# each SNP position refers to an anonymous hash of SNP alleles 
				# with their corresponding counts:
				$SNP_allele_freq_Aunat{$SNP_position}->{ $alleles_Aunat[$allele]->[$SNP_position] }++; 

			}
		}
##		# transform SNP allele counts into allele frequencies
		foreach my $SNP_position ( sort keys %SNP_allele_freq_Aunat ) {
			# print $columns[0], "/", $SNP_position+1, "\n";
			foreach my $allele ( keys %{ $SNP_allele_freq_Aunat{$SNP_position} } ) {
				$SNP_allele_freq_Aunat{$SNP_position}->{$allele} /= ($genotype_calls_Aunat*2);
			# 	print "\t", $allele, "\t", $SNP_allele_freq_Aunat{$SNP_position}->{$allele}, "\n"; 
			}
			print "\n";	
		} 

		
### Greixer		
		foreach my $ind_genotype (@Greixer_genotypes) {
			if ($ind_genotype ne "") { 
				$genotype_calls_Greixer++;
				# reformat homozygote genotypes from say AT to AT/AT
				if ($ind_genotype =~ /^([GATC]+)$/) {$ind_genotype = "$1/$1";} 
###				print $ind_genotype . "\n";
				# split genotype into it's two alleles
				my @alleles = split(/\//, $ind_genotype); 
###				print "@alleles",  "\n"; 
				foreach my $allele (@alleles) {
					push @alleles_Greixer, [ split(//, $allele) ];
				}
			}
		}
###			print $haplotypes_Aunat[0]->[0], "\n";
###			print @{$haplotypes_Aunat[0]}, "\n";
###			for my $hap (@haplotypes_Aunat) {
###				print "@{$hap}", "\n";
###			}

		# count SNP alleles for each SNP position	
		foreach my $SNP_position ( 0 .. $num_SNPs-1 ) {
			foreach my $allele ( 0 .. $#alleles_Greixer ) {
###				print $hap, " ", $SNP_position, "\n";
###				print $haplotypes_Aunat[$hap]->[$SNP_position], "\n";
				$SNP_allele_freq_Greixer{$SNP_position}->{ $alleles_Greixer[$allele]->[$SNP_position] }++; 
			}
		}
##		# transform SNP allele counts into allele frequencies
		foreach my $SNP_position ( sort keys %SNP_allele_freq_Greixer ) {
###			print $SNP_position+1, "\n";
			foreach my $allele ( keys %{ $SNP_allele_freq_Greixer{$SNP_position} } ) {
				$SNP_allele_freq_Greixer{$SNP_position}->{$allele} /= ($genotype_calls_Greixer*2);
###				print "\t", $allele, "\t", $SNP_allele_freq_Greixer{$SNP_position}->{$allele}, "\n"; 
			}
###			print "\n";	
		}
#		
#		
### now look for divergent SNPs between the two populations
### there are only bi-allelic SNPs in the export_sql_geno_13022012.tsv file
### the following code would not extract tri-allelic divergent SNPs! For this Nei's standard 
### genetic distance measure might be convenient
#		
		SNP: foreach my $SNP ( 0 .. $num_SNPs-1 ) {
				 # foreach allele that is found in the Greixer populations
				 foreach my $allele ( keys %{ $SNP_allele_freq_Greixer{$SNP} } ) {
					if ( $SNP_allele_freq_Greixer{$SNP}->{$allele} >= $threshold 
							&& ( !exists $SNP_allele_freq_Aunat{$SNP}->{$allele} # "!" with higher precedence than "not" is necessary here
								|| $SNP_allele_freq_Aunat{$SNP}->{$allele} <= 1-$threshold ) 
						) 

						
#							{ print "$columns[0]", "/", $SNP+1, "\n\t", "Greixer: ", 
#								$SNP_allele_freq_Greixer{$SNP}->{$allele}, "\n\t";
#								if (exists $SNP_allele_freq_Aunat{$SNP}->{$allele}) {
#									print "Aunat: ", $SNP_allele_freq_Aunat{$SNP}->{$allele}, "\n\n";}
#								else {print "Aunat: 0.0", "\n\n";}
#							next SNP;	
							#
							{ 
								#print STDERR "True\n";
								my $p_Greixer = $SNP_allele_freq_Greixer{$SNP}->{$allele};
								# initilize $p_Aunat		
								my $p_Aunat = 0;							  	

								if (exists $SNP_allele_freq_Aunat{$SNP}->{$allele}) {
									$p_Aunat = $SNP_allele_freq_Aunat{$SNP}->{$allele};
									print "True\n";
								}
							
								# initialize $F_st
								my $F_st = 0;
								my $p_mean = ($p_Greixer + $p_Aunat)/2;
								unless ( $p_mean == 1) { 
									$F_st = ( ($p_Greixer**2 + $p_Aunat**2)/2 - ($p_mean)**2 )
										/ ( $p_mean * (1-$p_mean) );
									}
									print "$columns[0]", "/", $SNP+1, "\n\t", "Greixer: ", $p_Greixer, "\n\t", 
									"Aunat: ", $p_Aunat, "\n\t", "Fst: ", $F_st, "\n\t", 
									"allele: ", $allele, "\n\n";
### Claudius: where is this formula taken from ?!
							
								print OUT join('	', @columns[0..11]), "\t", $genotype_calls_Aunat, 
								"\t", $genotype_calls_Greixer, "\t", $SNP+1, "\t", $F_st, "\t", 
								join('	', @columns[12..$#columns]), "\n";					
								next SNP; 
								# I am printing out the data for each SNP in the tag that meets the criteria
							}
							#elsif ( $SNP_allele_freq_Greixer{$SNP}->{$allele} <= 1-$threshold ) {
						#if ( exists $SNP_allele_freq_Aunat{$SNP}->{$allele} 
									#&& $SNP_allele_freq_Aunat{$SNP}->{$allele} >= $threshold )
								#{ 
										#my $p_Greixer = $SNP_allele_freq_Greixer{$SNP}->{$allele};
										#my $p_Aunat = $SNP_allele_freq_Aunat{$SNP}->{$allele};
										#my $p_mean = ($p_Greixer + $p_Aunat)/2;
										#my $F_st = ( ($p_Greixer**2 + $p_Aunat**2)/2 - ($p_mean)**2 )
										#/ ( $p_mean * (1-$p_mean) );
													#
													#print OUT join('	', @columns[0..11]), "\t", $genotype_calls_Aunat, "\t",
										#$genotype_calls_Greixer, "\t", $SNP+1, "\t", $F_st, "\t", 
												#join('	', @columns[12..$#columns]), "\n";
												#next SNP;
										#}
					 else {next SNP;}						
			} # close foreach allele
		} # close foreach SNP
} # close read line

system("tail -n +2 *_div_SNPs.tsv | sort -k 1n,1 -k 16nr,16  > temp") == 0 or die $!;
# system("cut -f 1 temp > cat_id") == 0 or die $!;
# system("cut -f 2- temp > rest") == 0 or die $!;
# system("paste rest cat_id > exchanged") == 0 or die $!;
# system("rm cat_id rest") == 0 or die $!;
# 
# open IN, "<", "exchanged" or die "Can't find file \"exchanged\".\n";
# $_ = <IN>;
# my @row = split(/\t/,$_);
# my $last_column = @row-1;
# close IN;
# 
# system("uniq -f $last_column < exchanged > sorted") == 0 or die $!; 
