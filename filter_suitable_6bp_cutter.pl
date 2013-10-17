#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: filter_suitable_6bp_cutter.pl
#
#        USAGE: ./filter_suitable_6bp_cutter.pl  
#
#  DESCRIPTION: This script searches a file containng 6bp cutters (retrieved from Rebase) for those RE's whose 
#				recognition site does not overlap the with SbfI restriction site. In a double-digest complexity reduction,
#				I don't want the two enzymes to interfere with each other. This script also filters out all enzymes with 
#				ambiguous recognition sequences.
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), claudiuskerth[at]gmail.com
# ORGANIZATION: Sheffield University
#      VERSION: 1.0
#      CREATED: 04/07/2012 15:33:20
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

my $usage = "
$0 <file_with_6bp_cutters_to_filter> 

This script searches a file containng 6bp cutters (e. g. retrieved from Rebase) for those RE's whose 
recognition site does not overlap the with SbfI restriction site. The input file should be in the format of
\"enzyme name\" \"any whitespace\" \"recognition sequence\" (without the quotation marks) with one line per enzyme. 
In a double-digest complexity reduction, I don't want the two enzymes to interfere with each other. 
This script also filters out all enzymes with ambiguous recognition sequences.
\n";
#The second input file is optional and should contain the 5bp barcodes that will be used with the RE's. 
#Even if ligation occurs after heat-inactivation of RE's, it might be better if the adapters don't contain
#recognition sequences for the RE's.

if (@ARGV == 0 || $ARGV[0] =~ /-h/) { die $usage; }

print "\n$0 ", "@ARGV", "\n";

my @bases = qw(A G C T);

my $SbfI = "CCTGCAGG";

# the following nested loops create all the combinations of 5bp sequences downstream (i. e. 3') 
# of the SbfI recognition sequence
my @seq = ();
for my $base (@bases) {
	my $nt_1 = $base;
	for my $base (@bases) {
		my $nt_2 = $base;
		for my $base (@bases) {
			my $nt_3 = $base;
			for my $base (@bases) {
				my $nt_4 = $base;
				for my $base (@bases) {
					my $nt_5 = $base;
					push @seq, join("", $nt_1, $nt_2, $nt_3, $nt_4, $nt_5, $SbfI, $nt_1, $nt_2, $nt_3, $nt_4, $nt_5);
				}
			}
		}
	}
}

# the input file should be a tab or space delimited list of enzyme names and their 
# recognition sequences
open(IN, "<", "$ARGV[0]") or die;

my %_6bp_cutters;

# store the enzymes names in a hash with recognition sequences as keys
# note: only one enzyme per recognition sequence will be kept
while (<IN>) {
	my ($enzyme, $recogn_seq) = $_ =~ /^(\w+)\s(\w+)/;
	# ignore lines with ambiguous recognition sequences
	next if $recogn_seq =~ /[^AGCT]/;
#	print $enzyme, " ", $recogn_seq, "\n";
	$_6bp_cutters{$recogn_seq} = $enzyme;
}
close IN;

# foreach possible sequence
foreach my $seq (@seq) {
	# foreach recognition sequence
	foreach my $recogn_seq (keys %_6bp_cutters) {
		# if the sequence contains the recognition sequence
		if ($seq =~ /$recogn_seq/) {
			print $_6bp_cutters{$recogn_seq}, "\tdoes overlap with SbfI recognition sequence.\n";
			# delete that enzyme from the hash
			delete $_6bp_cutters{$recogn_seq};
		}
	}
}
print "\n";

if ($ARGV[1]) {

# UNDER CONSTRUCTION
#	@seq = ();
#	# for checking the p2y adapter:
#	# this part should produce 64 sequences
#	my $p2y = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG";
#	
#	for my $base (@bases) {
#		my $nt_1 = $base;
#		for my $base (@bases) {
#			my $nt_2 = $base;
#			for my $base (@bases) {
#				my $nt_3 = $base;
#				for my $base (@bases) {
#					my $nt_4 = $base;
#					for my $base (@bases) {
#						my $nt_5 = $base;
#						push @seq, join("", $nt_1, $nt_2, $nt_3, $nt_4, $nt_5, $p2y);
#					}
#				}
#			}
#		}
#	}
	
	# check the P1 adapters for recognition sequences
	open(IN, "<", "$ARGV[1]") or die $!;
	my @barcodes = <IN>;
	chomp @barcodes;
	
	# illumina P1 adapter 5'-3'
	my $P1 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
	# reaminder of SbfI restriction site after adapter ligation
	my $sticky_end = "TGCAGG";
	my $seq = "";

	for my $barcode (@barcodes) {
		$seq = join("", $P1, $barcode, $sticky_end);
		# foreach recognition sequence
		foreach my $recogn_seq (keys %_6bp_cutters) {
			# if the sequence contains the recognition sequence
			if ($seq =~ /$recogn_seq/) {
				print $_6bp_cutters{$recogn_seq}, "\tcan cut in the P1 adapter with barcode $barcode.\n"; 
				# delete that enzyme from the hash
				delete $_6bp_cutters{$recogn_seq};
			}
		}
	}
}
print "\n";

# print out the reduce hash, now only containing enzymes whose recognition sequences can never overlap
# with the SbfI recognition sequence
print "6bp cutters with unambiguous recognition site,
no overlap with the SbfI restriction site 
and no cut in the P1 adapters are:", "\n";
foreach my $recogn_seq (sort { $_6bp_cutters{$a} cmp $_6bp_cutters{$b} } keys %_6bp_cutters) {
	print $_6bp_cutters{$recogn_seq}, "\t", $recogn_seq, "\n";
}
