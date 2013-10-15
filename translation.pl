#!/usr/bin/env perl
# translation.pl
use strict; use warnings;
use Bio::SeqIO;
use MyModule;

# 1. open input file given as exclusive command line argument
# 2. check fasta format
# 3. check IUPAC symbols
# 4. read in one fasta record at a time
# 5. send fasta sequence to MyModule::translation
# 6. modify MyModule::translation to search for longest ORF in all 6 reading frames
# 7. make this script report the length of the longest ORF, the number of the reading frame and the sequence header line

my $usage = "\n$0 <fasta file to translate>\n\n";
die $usage unless @ARGV;

my $infile = $ARGV[0];
my $seq_obj;

my $seqIO = Bio::SeqIO->new(-file => $infile,
							-format => "fasta"
							);
							

my $rf_1;
my $rf_2;
my $rf_3;
my $rf_4;
my $rf_5;
my $rf_6;
my @rf_1;
my @rf_2;
my @rf_3;
my @rf_4;
my @rf_5;
my @rf_6;
my %hash;
my $longest_ORF;

while($seq_obj = $seqIO->next_seq){
	print $seq_obj->display_id, "\n";
	print $seq_obj->desc, "\n";
	#print $seq_obj->translate->seq, "\n";
	#print $seq_obj->trunc(2, $seq_obj->length)->translate->seq, "\n";
	#print $seq_obj->trunc(3, $seq_obj->length)->translate->seq, "\n";
	#print $seq_obj->revcom->translate->seq, "\n";
	#print $seq_obj->revcom->trunc(2, $seq_obj->length)->translate->seq, "\n";
	#print $seq_obj->revcom->trunc(3, $seq_obj->length)->translate->seq, "\n";
	$rf_1 = $seq_obj->translate(-frame => 0)->seq;
	@rf_1 = split('\*', $rf_1);
	$hash{ longest_element(@rf_1) } = "reading frame 1";
	#foreach(@rf_1){ print length($_), "\n";}
	$rf_2 = $seq_obj->translate(-frame => 1)->seq;
	@rf_2 = split('\*', $rf_2);
	$hash{ longest_element(@rf_2) } = "reading frame 2";
	$rf_3 = $seq_obj->translate(-frame => 2)->seq;
	@rf_3 = split('\*', $rf_3);
	$hash{ longest_element(@rf_3) } = "reading frame 3";
	$rf_4 = $seq_obj->revcom->translate(-frame => 0)->seq;
	@rf_4 = split('\*', $rf_4);
	$hash{ longest_element(@rf_4) } = "reading frame 4";
	$rf_5 = $seq_obj->revcom->translate(-frame => 1)->seq;
	@rf_5 = split('\*', $rf_5);
	$hash{ longest_element(@rf_5) } = "reading frame 5";
	$rf_6 = $seq_obj->revcom->translate(-frame => 2)->seq;
	@rf_6 = split('\*', $rf_6);
	$hash{ longest_element(@rf_6) } = "reading frame 6";

	$longest_ORF = longest_element( keys(%hash) );
	print $hash{ $longest_ORF }, "\n";
	print $longest_ORF, "\t", $seq_obj->display_id, "\t", $seq_obj->desc, "\n";
	print length( $longest_ORF ) ;
	print "\n\n";
}
#####################################
# SUBROUTINE
# name: longest_element
# receives: array of outative ORF's
# returns: longest ORF in array
# http://stackoverflow.com/questions/4182010/the-fastest-way-execution-time-to-find-the-longest-element-in-an-list
#####################################
sub longest_element{
	my $max = -1;
	my $max_ref;
	for(@_){
		if(length($_) > $max){
			$max = length($_);
			$max_ref = \$_;
		}
	}
	return $$max_ref;
}
#####################################

#
#my %fasta_record;
#
#while(<>){
#	chomp;
#	if(/^>/){ $header = $_; }
#	else{ $seq .= $_; }
#
#}
#
##open(my $in, "<", $seq) || die("can't open datafile: $!");
#
#
#my $IUPAC = "acgtwsmkrybdhvn";
#die "non-IUPAC symbols used" if ($seq =~ m/[^$IUPAC]/i);
#
##die "non-DNA character detected. Only use A, G, C or T, please\n" if ($seq =~ /[^AGCT]/i);
#
##if ($seq =~ /agc/i) {$seq .= "this is a test";}
##print "$seq\n";
# 
#my @protein = MyModule::translation($seq);
#for (my $i = 0; $i < @protein; $i++) {
#	print $i+1, ". reading frame: ", "$protein[$i]\n";
#}
#
