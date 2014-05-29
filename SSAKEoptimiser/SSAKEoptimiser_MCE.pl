#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: SSAKEoptimiser_MCE.pl
#
#        USAGE: ./SSAKEoptimiser_MCE.pl  
#
#  DESCRIPTION: This script runs the SSAKE assembler on an array of input files
#  				for a range of kmer lengths and keeps the output containing the
#  				the longest contig for each input file. It does the SSAKE runs 
#  				in parallel with the Perl module Many-core Engine (MCE), maximising 
#  				the usage of available cores.
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 15/05/14 12:12:18
#     REVISION: ---
#===============================================================================


use strict;
use warnings;

use Time::HiRes qw(time sleep);

my $start = time();

my $usage = "
$0 -in *fa

-in     takes input file names (shell character expansion allowed)
-spare  specify the number of cores that this script should NOT use
-h      prints this help
\n";

#-------------------------------------------------------------------------------
#  get input and sort by size
#-------------------------------------------------------------------------------

my @in_filenames;
my $spare = 0;

parse_command_line(@ARGV);

die $usage unless @in_filenames;
for(@in_filenames){
	unless (-e $_){ die "file does not exist: $_\n"}
}


use MCE::Loop chunk_size => 1;

MCE::Loop::init {
	max_workers => 10
};

my %h;
# get file size for each input file
%h = mce_loop {
	my $filename = $_;
	my $filesize = -s $filename;
	MCE->gather($filename, $filesize);
} @in_filenames;

## sort files by size descendingly ##
my @in_filenames_sorted = sort {$h{$b} <=> $h{$a}} keys %h;


#-------------------------------------------------------------------------------
#  run SSAKE in parallel
#-------------------------------------------------------------------------------

use MCE::Step;
use MCE::Util 'get_ncpu';

my $n_cores = get_ncpu();
my $spawn = $n_cores- 2 - $spare;
die "No cores left for SSAKE run\n" if ($spawn <= 0);
print "Using ", $spawn+2, " of $n_cores cores on this machine.\n";

## create reference to input
my $input =  [@in_filenames_sorted[30 .. 40]]; 

## run MCE Step ##
mce_step {
	task_name => ['command', 'SSAKE'],
	max_workers => [ 1, $spawn ], # main: 1, command: 1, SSAKE: ncpu-2-spare
	chunk_size => 1,
	gather => gather_output()
}, \&setup_commandline, \&run_SSAKE, $input;

#printf STDERR "\n## Comnpute time: %0.3f secs\n\n", time() - $start;
printf STDOUT "%d\t%.3f\n", $spawn+2, time() - $start;


#-------------------------------------------------------------------------------
#  subroutines
#-------------------------------------------------------------------------------

sub parse_command_line{
	while($_ = shift){
#		print "$_\n";
		if(/^-+in/){
			while($_ = shift) {
				last if /^-+/;
				push @in_filenames, $_;
			}
		}
		if(defined($_)){
			if(/^-+h/){ die $usage }
			elsif(/^-+spare/){ $spare = shift; }
		}
	}
}

## This function is run by the manager process ##
sub gather_output {
	# this hash temporarily stores file names and kmer lengths
	# of completed SSAKE runs
	my %tmp; 
	return sub {
		my ($file) = @_;
		$tmp{$file}++;
		if($tmp{$file} == 33-11+1){
#			MCE->print("Finding longest contig for $file.\n");
			find_longest_contig($file);
			delete $tmp{$file};
		}
		return;
	}
}

sub setup_commandline {
	my $file = $_;
	for my $kmer (11..33){
		MCE->step($file, $kmer);
	}
}

sub run_SSAKE {
	my ($mce, $file, $kmer) = @_;
	my $cmd = "SSAKE -f $file -w 1 -o 1 -m $kmer -c 1 > /dev/null";
#	MCE->print("Assembling $file with kmer length of ", $kmer, "\n");
	# this does a fork for the command,
    # while the worker waits for the forked process to complete
	system($cmd) == 0 or die $!; 
	MCE->gather($file);
}

#===  FUNCTION  ================================================================
#         NAME: Find_Longest_Contig
#      PURPOSE: Finding the output file containing the longest contig among all
#      			output files produced by different SSAKE runs with different
#      			kmer sizes. Removes all other output files (clean-up).
#   PARAMETERS: input file name
#      RETURNS: nothing (I think)
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub find_longest_contig {
	my $file = shift;

	# get all output file names containing the contigs
	my @output = `ls $file*contigs`; 
	chomp @output;

	my %longest_contig_hash;

	# get the length of the longest contig for each SSAKE x kmer run
	for my $output (@output){

		if(-z $output){ # if the contig file is empty
			print("Contig file is empty for $output\n");
			$output =~ s/contigs$/\*/; 
			system("rm -f $output")==0 or die $!;
			next; 
		}

		my $longest_contig_length = 0; # reset to 0 for each kmer output file
		open my $IN, "<", $output or die $!;
		my $length;
		while(<$IN>){
			next if /^>/;
			chomp;
			$length = length($_);
			$longest_contig_length = $length if ($length > $longest_contig_length);
		}
		close $IN;

#		print("$output : $longest_contig_length\n");
		unless ($longest_contig_length){
			print("Found no contig length for $output\n");
			next;
		}; 
		$longest_contig_hash{$output} = $longest_contig_length;
	}

	# exit if no contigs could be assembled
	unless(keys %longest_contig_hash){
		print("No contigs could be assembled for $file.\n");
		map { s/contigs$/\*/ } @output;
		system("rm -f @output")==0 or die $!;
		return;
	}

	# sort output file names by contig length
	my @remove = sort {$longest_contig_hash{$b} <=> $longest_contig_hash{$a}} keys %longest_contig_hash;

	# save the name of the output file containing the longest contig from
	# all SSAKE runs
	my $keep = shift(@remove);
#	print("File containing longest contig:\n", "$keep : $longest_contig_hash{$keep}\n\n");

	map { s/contigs$/\*/ } @remove;
	system("rm -f @remove")==0 or die $!;

	return;
}

