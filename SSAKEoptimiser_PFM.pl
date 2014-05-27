#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: SSAKEoptimiser_PFM.pl
#
#        USAGE: ./SSAKEoptimiser_PFM.pl  
#
#  DESCRIPTION: This script uses Parallel:ForkManager to parallise de novo
#  				assemblies with SSAKE.
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 23/05/14 18:40:10
#     REVISION: ---
#===============================================================================

 # :TODO      :27/05/14 13:08:18:CEK: check for required packages
 # :TODO      :27/05/14 13:08:51:CEK: parse command line arguments for max forks and input file path

use strict;
use warnings;

use Time::HiRes qw(time sleep);

my $start = time();

#-------------------------------------------------------------------------------
#  get input and sort by size
#-------------------------------------------------------------------------------

## get input file names
my @files = glob("input/*.fa");

## get line count for each file
my %h;
for (@files){
	my $filename = $_;
	my $filesize = -s $filename;
	$h{$filename} = $filesize;
}
## sort files by line count descendingly
my @files_sorted = sort {$h{$b} <=> $h{$a}} keys %h;

use Parallel::ForkManager;

# initialise PFM object with maximum number of processes to fork
my $pm = new Parallel::ForkManager(22);

#-------------------------------------------------------------------------------
#  gather function
#-------------------------------------------------------------------------------
my %gather;

$pm->run_on_finish(
	## this code block is run in the parent process,
	## when a child has terminated (finished).
	sub { 
		## the third parameter received by the 'run_on_finish' method
		## is the 'identification of the process'
		## given when starting the child process
		## see Callbacks in the documentation of Parallel::ForkManager
		my $file = $_[2];
		$gather{$file}++;
		if($gather{$file}==33-11+1){
			print("Finding longest contig for $file.\n");
			find_longest_contig($file);
			delete $gather{$file};
		}
	}
);

#-------------------------------------------------------------------------------
#  run SSAKE in parallel
#-------------------------------------------------------------------------------

## run parallel assemblies
for my $file (@files_sorted[30..40]){
	for my $kmer (11..33){
		print("Assembling $file with kmer length of ", $kmer, "\n");
		my $pid = $pm->start($file) and next;
			my $cmd = "SSAKE -f $file -w 1 -o 1 -m $kmer -c 1 > /dev/null";
			system($cmd) == 0  or die $!;
			$pm->finish;
	}
}
$pm->wait_all_children;

printf STDERR "\n## Comnpute time: %0.3f secs\n\n", time() - $start;

#-------------------------------------------------------------------------------
#  subroutines
#-------------------------------------------------------------------------------

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

	my $longest_contig_length = 0;
	my %longest_contig_hash;

	# get the length of the longest contig for each SSAKE x kmer run
	for my $output (@output){

		if(-z $output){ # if the contig file is empty
			print("Contig file is empty for $output\n");
			$output =~ s/contigs$/\*/; 
			system("rm -f $output")==0 or die $!;
			next; 
		}

		open my $IN, "<", $output or die $!;
		my $length;
		while(<$IN>){
			next if /^>/;
			chomp;
			$length = length($_);
			$longest_contig_length = $length if ($length > $longest_contig_length);
		}
		close $IN;
		print("$output : $longest_contig_length\n");
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

	my @remove;
	# sort output file names by contig length
	foreach my $output (sort {$longest_contig_hash{$b} <=> $longest_contig_hash{$a}} keys %longest_contig_hash){
		push(@remove, $output);
	}


	# save the name of the output file containing the longest contig from
	# all SSAKE runs
	my $keep = shift(@remove);
	print("File containing longest contig:\n", "$keep : $longest_contig_hash{$keep}\n\n");


	map { s/contigs$/\*/ } @remove;
	system("rm -f @remove")==0 or die $!;

	return;
}

