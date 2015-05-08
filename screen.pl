#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: screen.pl
#
#        USAGE: ./screen.pl  
#
#  DESCRIPTION: This script takes sam records of SE reads from read pairs where
#               both reads mapped but to different reference contigs, i. e. discordantly:
#               samtools view -f64 -F2 file.sam
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 25/03/15 16:43:28
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

my (@fields, $name, %contigs, $length, $flag, $SE_contig, $PE_contig, $SE_mappos, $PE_mappos, $contig_length, $SE_OK, $PE_OK);

my $SE_read_length = 46;
my $PE_read_length = 51;
#my $SE_read_length = 70;
#my $PE_read_length = 70;
#my $SE_read_length = 0;
#my $PE_read_length = 0;

while(<>){
	#----------------------------------------------
	# read in contig lengths from bam header
	#----------------------------------------------
	if(/^\@/){ # if this is a header line
		if(/^\@SQ/){
			@fields = split;
			($name = $fields[1]) =~ s/^SN://;
			($length = $fields[2]) =~ s/^LN://;
			$contigs{ $name } = $length;
		}
		next;
	}

	@fields = split;
	($flag, $SE_contig, $SE_mappos, $PE_contig, $PE_mappos) = @fields[1,2,3,6,7];
#	$SE_read_length = length($fields[9]);
#	$PE_read_length = 96 - $SE_read_length + 4;
#	print join("	", ($flag, $SE_contig, $SE_mappos, $PE_contig, $PE_mappos)), "\n";
   #--------------------------------------------------
   # check if the read pair would have had enough 
   # space on either reference contig to map properly
   #--------------------------------------------------
	($SE_OK, $PE_OK) = (0, 0);
	# if SE read is mapped as revcomp
	if($flag & 16){
		# is there enough space upstream for any possible paired end read properly mapping on the same contig?
		if( ($SE_mappos + $SE_read_length -1) >= 900 ){ $SE_OK = 1 }; # 900 is the maximum possible fragment size I gave to bowtie2
	# if the SE read is mapped on same strand as reference contig sequence
	}else{
		die "could not find length of contig $SE_contig\n" unless exists $contigs{ $SE_contig };
		$contig_length = $contigs{ $SE_contig };
		# is there enough space downstream for any properly mapping paired end read?
		if( ($contig_length - $SE_mappos +1) >= 900 ){ $SE_OK = 1 }; 	
	}
	# if PE read is mapped as revcomp
	if($flag & 32){
		if( ($PE_mappos + $PE_read_length -1) >= 900 ){ $PE_OK = 1 };
	# if the PE read is mapped on same strand as reference contig sequence
	}else{
		die "could not find length of contig $PE_contig\n" unless exists $contigs{ $PE_contig };
		$contig_length = $contigs{ $PE_contig };
		if( ($contig_length - $PE_mappos +1) >= 900 ){ $PE_OK = 1 }; 	
	}

	#-------------------------------------------------------------
	# print out the sam record if either check turned out ok
	# indicating that at least one of the two reference contigs
	# could have provided the reference for both reads in the pair
	#-------------------------------------------------------------
	if( $SE_OK || $PE_OK ){
		print join("	", @fields), "\n";
	}
}

#foreach (keys %contigs){
#	print $_, "\t", $contigs{ $_ }, "\n";
#}

#if(/^\@SQ/){@l = split; $l[1] =~ s/^SN://; $l[2] =~ s/^LN://; $H{$l[1]} = $l[2];} END{foreach (keys %H){ print $_, "\t", $H{$_}, "\n";}}
