#!/usr/bin/env perl

use warnings;
use strict;

use Config;
$Config{useithreads} or die "Recompile Perl with threads to run this program.";

use threads qw(stringify);
use Thread::Semaphore;


my $usage="
$0 

Options:

-c   number of cores to use (default 1)

This is a wrapper script for SSAKE.
It takes all the files ending with \'fa\' in the current directory
and runs SSAKE for a range of kmer values from 11 to 33.
It keeps the assembly with the longest contig for each input file.
It takes the number of cores to use as its single argument.
\n";

my $ncores = 1;

parse_command_line();

print "Using up to $ncores cores in parallel.\n";

my $sem = Thread::Semaphore->new($ncores);

# get input file names
my @infiles = glob("*fa");
chomp @infiles;

my @main_threads;
# spawn a new thread for each input file
foreach my $infile (@infiles){
	my $thr = threads->new(\&optimise_SSAKE, $infile);
	push @main_threads, $thr;
}

$_->join for @main_threads; # wait for everything to finish


########## SUBROUTINES #############

sub optimise_SSAKE {

	$sem->down;
	my $input = shift;
	my @threads;
	my @returned;
	my @existing_output;

	# run SSAKE on input file for a range of kmer lengths
	KMER: for(my $kmer_length=11; $kmer_length <= 33; $kmer_length++){
		while(<$input*_m$kmer_length*contigs>){
			unless (-z $_){ print  "\nSkipping SSAKE run on $input for kmer length $kmer_length. Output already exists.\n"; next KMER; }
		}
#		$threads{$kmer_length} = threads->new( 
		$sem->down;
		my $thr = threads->new( 
		sub	{ 
				print  "Running SSAKE on $input with kmer length $kmer_length.\n";
				my ($input, $kmer_length) = @_; 
				if (system("SSAKE -f $input -w 1 -o 1 -m $kmer_length -c 1 > /dev/null")){$sem->up; return $!;};
				$sem->up;
				return "Finished assembling reads for $input with kmer length $kmer_length.";
			},
			($input, $kmer_length)
		);
#		print "Thread $thr started ...\n";
		push @threads, $thr;
	}

	print $_->join, "\n" for (@threads);

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
		return "No contigs could be assembled for $input.\n";
	}

	my @remove;
	# sort output file names by contig length
	foreach my $output (sort {$longest_contig_hash{$b} <=> $longest_contig_hash{$a}} keys %longest_contig_hash){
		push(@remove, $output);
	}


	# remove the name of the output file containing the longest contig from
	# all SSAKE runs
	$_ = shift(@remove);
	print "File containing longest contig:\n", "$_\n\n";


	# remove all other SSAKE output files
	my $file;
	while(@remove){
		$file  = shift @remove;
		$file =~ s/contigs$//;
		system("rm -f $file*")==0 or die $!;
	}
	$sem->up;
}


sub parse_command_line {
	while(@ARGV){
		$_ = shift @ARGV;
		if(/-c/){ $ncores = shift @ARGV; }
		elsif(/-h|--h/){ die $usage; }
		else{ die $usage }
	}
}
