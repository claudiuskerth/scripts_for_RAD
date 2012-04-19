#!/usr/bin/env perl
use strict; use warnings;

my $usage = "\n$0 <file names to rename> <barcode_ind file>

<file names to rename> can include wildcard expression which needs to be in double quotes
<barcode_ind file> should be a tab delimited file in the format of:
			AGCTA	30-1
			ATCGT	30-2
			etc.
-n no action, just print what would be done 
\n";

die $usage unless @ARGV >= 2;

my @files = glob($ARGV[0]) or die $!;

#foreach my $file (@files) { print $file, "\n"; }

open(IN, "<", $ARGV[1]) or die $!;

my ($barcode, $ind, %barcode_ind, $no_action);

if ( defined($ARGV[2]) && $ARGV[2] =~ /-n/) { $no_action = "true"; }

while(<IN>){
	chomp;
	($barcode, $ind) = split /\t/, $_;
	next unless ( defined($barcode) && defined($ind) );
	$barcode_ind{$barcode} = $ind;
}

#foreach my $barcode (sort keys %barcode_ind) {
#print $barcode, "\t", $ind, "\n";
#}

File: 
foreach my $file (@files) {
	foreach my $barcode ( sort keys %barcode_ind ) {
		if ($file =~ /$barcode/) {
			(my $file_new = $file) =~ s/$barcode/$barcode_ind{$barcode}/;
			#print $file, "\t", $file_new, "\n";
			if ($no_action) { print "$file will be renamed as $file_new\n"; }
			else {
				system("mv $file $file_new") == 0 or die "cannot rename\n";
			}
			next File;
		}
	}
}
