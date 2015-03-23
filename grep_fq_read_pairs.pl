#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: collect_SbfI_reads.pl
#
#        USAGE: grep_fq_read_pairs.pl
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
#      CREATED: 05/03/15 13:50:30
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

my ($SE_file, $PE_file, $search_seq);

$search_seq = shift @ARGV;

# provide SE file names as arguments
foreach my $filename (@ARGV){
	# open SE file with filename, unzip if necessary
	if($filename =~ /gz$/){
		open( $SE_file, "zcat $filename |") or die $!;
	}else{
		open( $SE_file, $filename ) or die $!;
	}
	# open corresponding PE file
	$filename =~ s/fq_1/fq_2/;
	if($filename =~ /gz$/){
		open( $PE_file, "zcat $filename |") or die $!;
	}else{
		open( $PE_file, $filename ) or die $!;
	}

	# open output files
	$filename =~ s/(.*)(\.fq_2).*/$1_SbfI_reads$2/;
	open(PE_OUT, ">", $filename) or die $!;
	$filename =~ s/fq_2/fq_1/;
	open(SE_OUT, ">", $filename) or die $!;

	while( my @s = fq_record($SE_file), my @p = fq_record($PE_file) ){
		if($s[1] =~ /$search_seq/ or $p[1] =~ /$search_seq/){
			print SE_OUT @s;
			print PE_OUT @p;
		}
	}
	close SE_OUT; 
	close PE_OUT;
}

sub fq_record {
	my $fh = shift;
	my $head = <$fh>;
	my $seq = <$fh>;
	my $qh = <$fh>;
	my $qual = <$fh>;
	if( $head && $seq && $qh && $qual ){
		return($head, $seq, $qh, $qual);
	}
	else{ return } # note this return without argument is necessary, otherwise $qual gets returned automatically
}
