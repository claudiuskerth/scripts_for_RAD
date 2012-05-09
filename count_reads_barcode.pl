#!/usr/bin/env perl
use strict; use warnings;
#
#==========================================================================
#  This program is free software; you can redistribute it and/or modify 
#  it under the terms of the GNU General Public License as published by 
#  the Free Software Foundation; either version 2 of the License, or    
#  (at your option) any later version.                                  
#==========================================================================
#
#open(IN, "<", $ARGV[0]) or die $!;

# pipe the read input into this programme
# specify file containing used barcodes sequences with "-b" switch
#

die "no command line arguments\n" unless @ARGV == 2;

my $bc;
while(@ARGV) {
	$_ = shift @ARGV;
	if ( $_ =~ /-b/ ) { $bc = shift @ARGV; }
}

open(BC, "<", $bc )or die $!;

my @barcodes = <BC>;
close BC;
foreach my $barcode (0..$#barcodes) { chomp($barcodes[$barcode]); } 
#foreach my $barcode (@barcodes) { print $barcode, "\n"; }

# check that barcodes have at least two mismatches
# between each other
my ($mismatch, $num_comp);
for ( my $i = 0; $i < @barcodes-1; $i++ ) {
	for (my $j = $i+1; $j < @barcodes; $j++ ) {
		$mismatch = ($barcodes[$i] ^ $barcodes[$j]) =~ tr[\001-\255][];
		#print "$mismatch\n";
	       	if ( $mismatch < 2 ) { 
			die "not enough mismatch between barcodes!\n", 
			$barcodes[$i], " ", $barcodes[$j], "\n";
		}
	}
}

my %barcode_count;
my $seq = "";
my $diff = 0;
my $no_barcode = 0;
my $restr_site = "";
my $no_radtag = 0;
my $sum = 0;
my $read_count = 0;


Line:
while(<>){
	next if ( ($.-2) % 4 );
	$read_count++;
#	print "$_\n";
	$seq = substr $_, 0, 5;
#	print $seq, "\n";
	foreach my $barcode (@barcodes) {
		$diff = ($seq ^ $barcode) =~ tr[\001-\255][];
#		print "$diff\n";
		# allow for one mismatch to the given barcode sequence:
		if ( $diff <= 1 ) { 
			$barcode_count{$barcode}->{total}++; 
			#next Line;
			$restr_site = substr($_, 5, 6);
			# allow for one mismatch to the remainder of the SbfI sequence:
			$diff = ($restr_site ^ "TGCAGG") =~ tr[\001-\255][];
			if ( $diff <= 1 ) {
				next Line;
			}
			$barcode_count{$barcode}->{no_radtag}++;
			print STDERR "$restr_site\n";
			$no_radtag++;	
			next Line;
		}
	}
	$no_barcode++;
}

# report
print "barcode\ttotal\tno radtag\n";

foreach my $barcode ( sort keys %barcode_count ) {
	print $barcode, "\t", $barcode_count{$barcode}->{total}, "\t", $barcode_count{$barcode}->{no_radtag}, "\n";
	$sum += ($barcode_count{$barcode}->{total} - $barcode_count{$barcode}->{no_radtag});
}
print "total with barcode and remainder of restriction site: ", $sum, "\n";
print "read ", $read_count, " sequences\n";
print "no barcode: ", $no_barcode, "\n";
print "with barcode but no rad tag: ", $no_radtag, "\n";
#close IN;

