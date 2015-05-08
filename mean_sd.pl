#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: mean_sd.pl
#
#        USAGE: ./mean_sd.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 07/05/15 13:53:07
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

my ($c, $ssq, $mean, $var, $sd);

while(<>){
	next if /^$/;
	chomp;
	die "found non-digit character in input\n" if /\D/;
	$c += $_;
	$ssq += $_**2;
}
$mean = $c/$.;
$var = $ssq/$. - $mean**2;
$sd = $var**(1/2);
printf "%.1f\t%.2f\n", $mean,  $sd; 
