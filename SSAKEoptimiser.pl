#!/usr/bin/env perl
use warnings;
use strict;

my $usage="
$0 <input file for SSAKE>

This is a wrapper script for SSAKE.
It iterates through an 11..33 kmer range and keeps the
assembly with the longest contig. Use GNU parallel
if you want to run this script on many input files in parallel.
\n";

die $usage unless @ARGV;

# get input file names
my $input = $ARGV[0];
chomp $input;

# run SSAKE on input file for a range of kmer lengths
for(my $kmer_length=11; $kmer_length <= 33; $kmer_length+=1){
	system("SSAKE -f $input -w 1 -o 1 -m $kmer_length -c 1 > /dev/null")==0 or die $!;
}

# get output file names containing the contigs
my @output = `ls $input*contigs`;

my $longest_contig_length = 0;
my %longest_contig_hash;

# get the length of the longest contig for each SSAKE run
for my $output (@output){
	chomp $output;
	if(-z $output){ # if the contig file is empty
		$output =~ s/contigs$//; 
		system("rm -f $output*")==0 or die $!;
		next; 
	}
	my $cmd = "awk \'NR%2==0\' " . $output . " | perl -ne\'chomp; print length, \"\\n\";\' | sort -rn | head -n1";
	die "$output\n" unless ($longest_contig_length = `$cmd`);
	die "Found no contig length\n" unless ($longest_contig_length); 
	chomp $longest_contig_length;
	$longest_contig_hash{$output} = $longest_contig_length;
}

# exit if no contigs could be assembled
unless(keys %longest_contig_hash){
	print "No contigs could be assembled for $input.\n";
	exit;
}

my @remove;
# sort output files names by contig length
foreach my $output (sort {$longest_contig_hash{$b} <=> $longest_contig_hash{$a}} keys %longest_contig_hash){
	print $output, "\t", $longest_contig_hash{$output}, "\n";
	push(@remove, $output);
}

print "\n";

# remove the name of the output file containing the longest contig from
# all SSAKE runs
$_ = shift(@remove);
print "File containing longest contig:\n", "$_\n\n";


# remove all other SSAKE output files
print "Removing files:\n";
while(@remove){
	my $file  = shift @remove;
	$file =~ s/contigs$//;
    print "$file*\n";
	system("rm -f $file*")==0 or die $!;
}

