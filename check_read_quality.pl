#!/usr/bin/env perl
use strict; use warnings;

# To do:
# - make script aware of different quality score encodings
# - make the script read in .fastq format data, currently
#   it expects each line to contain quality scores


my $usage = "
This programme checks whether input lines of quality scores
fall below a given threshold by sliding a window down each line
and calculating the average quality score at each step.
It reports when it finds a line in which the average quality score
in the window fell below the given threshold. 

$0

-w sliding window_size, must be an integer
-t quality threshold for the mean quality in the window
-enc [33] for +33 encoding, [64] for +64 encoding of qualtiy scores
\n";

my ( $window, $threshold, $qual, $encod );

parse_command_line();
print STDERR $window, "\t", $threshold, "\n";

die "the window size must be a positive interger above 0\n" 
	if ( ($window < 1) or ($window - int($window) > 0) );
die "the quality threshold must be between 0 and 60\n"
	if ( $threshold < 0 or $threshold > 60 );

# read in from a stream of quality score lines	
while(<>){
	chomp;
	$qual = 0;
	# first window:
	for (my $i = 0; $i < $window; $i++){
		$qual += ord( substr($_, $i, 1) ) - $encod;
	}

	if ( $qual < $threshold*$window ){ 
	# Don't divide $qual by $window and compare it to $threshold.
	# This would lead to inaccurate results.
		printf "%s%s%s%.19f%s", "found read with below threshold quality:\n",
			substr($_, 0, $window), "\n",
			$qual/$window, "\n";
	}
	# sliding window:
	for (my $i = $window; $i < length($_); $i++){
		$qual -= ( ord( substr($_, $i-$window, 1) ) - 33 );
		$qual += ( ord( substr($_, $i, 1) ) - 33 );
		
		if ( $qual < $threshold*$window ){
			printf "%s%s%s%.19f%s", "found read with below threshold quality:\n",
				substr($_, $i-$window+1, $window), "\n",
				$qual/$window, "\n";
		}	
	}
}

sub parse_command_line{
	die $usage unless @ARGV;
	while (@ARGV){
		$_ = shift @ARGV;
		if ( $_ =~ /^-w$/ ) { $window = shift @ARGV; }
		if ( $_ =~ /^-t$/ ) { $threshold = shift @ARGV; }
		if ( $_ =~ /^-enc$/ ) { $encod = shift @ARGV; }
	}

}
