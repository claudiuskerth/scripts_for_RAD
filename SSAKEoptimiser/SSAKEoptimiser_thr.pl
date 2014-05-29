#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: SSAKEoptimiser_thr.pl
#
#        USAGE: ./SSAKEoptimiser_thr.pl  
#
#  DESCRIPTION: This script opimises the SSAKE de novo assembler in parallel
#  				using Perls ithreads.
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 27/05/14 13:27:05
#     REVISION: ---
#===============================================================================

use strict;
use warnings;


use Time::HiRes qw(time sleep);

my $start = time();

my $usage = "
$0 -in *fa

-in     takes input file names (shell character expansion allowed)
-spare  specify the number of cores that this script should NOT use (default 0)
-h      prints this help
\n";


#-------------------------------------------------------------------------------
#  parse command line 
#-------------------------------------------------------------------------------

my @in_filenames;
my $spare = 0;

parse_command_line(@ARGV);

die $usage unless @in_filenames;
for(@in_filenames){
	unless (-e $_){ die "file does not exist: $_\n"}
}

#-------------------------------------------------------------------------------
#  get input and sort by size
#-------------------------------------------------------------------------------

### get input file names
#my @files = glob("input/*.fa");

## get line count for each file
my %h;
for (@in_filenames){
	my $filename = $_;
	my $filesize = -s $filename;
	$h{$filename} = $filesize;
}
## sort files by line count descendingly
my @in_filenames_sorted = sort {$h{$b} <=> $h{$a}} keys %h;



#-------------------------------------------------------------------------------
#  set number of cores to use
#-------------------------------------------------------------------------------
use MCE::Util 'get_ncpu';

my $n_cores = get_ncpu();
my $max_threads = $n_cores - 1 - $spare; # the manager thread spawns worker threads
die "No cores left for SSAKE run\n" if ($max_threads <= 0);
print "Using ", $max_threads+1, " of $n_cores cores on this machine.\n";

#-------------------------------------------------------------------------------
#  run SSAKE
#-------------------------------------------------------------------------------
use threads 'stringify';
use threads::shared;
use Thread;
use Thread::Semaphore;

my $mutex = Thread::Semaphore->new($max_threads);

my %gather :shared;

for my $file (@in_filenames_sorted[30..40]){
	for my $kmer (11..33){
		$mutex->down;
#		print("Assembling $file with kmer length of ", $kmer, "\n");
		my $thr = threads->create('run_SSAKE', $file, $kmer);
#		print "Thread $thr started ...\n";
	}
}

$_->join for threads->list();

#printf STDERR "\n## Comnpute time: %0.3f secs\n\n", time() - $start;
printf STDOUT "%d\t%.3f\n", $max_threads+1, time() - $start;

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

sub run_SSAKE {
	my ($file, $kmer) = @_;
	my $cmd = "SSAKE -f $file -w 1 -o 1 -m $kmer -c 1 > /dev/null";
	system($cmd)==0 or die $!;
		{
			lock %gather;
			gather($file);
		}
	$mutex->up;
}

sub gather {
	my ($file) = @_;	
	$gather{$file}++;
	if($gather{$file} == 33-11+1){
#		print("Finding longest contig for $file.\n");
		find_longest_contig($file);
		delete $gather{$file};
	}
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
