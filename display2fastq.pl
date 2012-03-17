#!/usr/bin/perl
use strict; use warnings;
use Getopt::Std;
our ($opt_l, $opt_v, $opt_h);
getopts('vhl:');

my $version = "1.0, 15/08/2011, author: Claudius Kerth\n";

my $usage = "
This script displays reads in two fastq files and displays headers or sequences of each file side by side,
i. e. first header from file 1 left of first header of file 2 and so on. With the option <check> it will
check whether the headers of both input files match, which can be used to check whether a single end fastq
file matches a paired-end fastq file.
usage: display2fastq.pl [options] <arguments...>
options:
	-v version
	-h help
	-l <header> or <seq>, to display header or sequence; <check> to check if headers match
	
";

if ($opt_v) { print $version; exit;
}

if ($opt_h) { print $usage; exit;}


my ($single_file, $paired_file) = @ARGV;

open my $single, "<", $single_file or die;
open my $paired, "<", $paired_file or die;

if (!$opt_l){ print $usage; exit;}
elsif ($opt_l eq 'seq') {
# print out single and paired end reads next to each other
	while (my $s = <$single>, my $p = <$paired>){
		if ($. == 2) {
			chomp $s;
			print "$s $p";
		}
		elsif ( ($. - 2) % 4 == 0 ) {
			chomp $s;
			print "$s $p";
		}
    	}
exit;
} 
elsif ($opt_l eq 'header') {
# print headers of single and paired end reads next to each other
	while (my $s = <$single>, my $p = <$paired>){
		if ($. == 1) {
		   chomp $s;
		   print "$s\t$p";
		}
		elsif ( ($. - 1) % 4 == 0 ) {
			chomp $s;
			print "$s\t$p";
		}
	}
exit;
}
elsif ($opt_l eq 'check') {
# print headers of single and paired end reads next to each other
	while (my $s = <$single>, my $p = <$paired>){
		if ($. == 1) {
			$s = substr($s,0,-2);
		   	$p = substr($p,0,-2);
#		    	print "$s\t$p\n";
			unless ( $s eq $p ) { print "non-matching headers on input line $.\n"};
			 
		}
		elsif ( ($. - 1) % 4 == 0 ) {
		   	$s = substr($s,0,-2);
		   	$p = substr($p,0,-2);
#		   	print "$s\t$p\n";
			unless ( $s eq $p ) { print "non-matching headers on input line $.\n"};
		}
	}
exit;
}

close $single;
close $paired;

