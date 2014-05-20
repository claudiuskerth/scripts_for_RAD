#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: SSAKEoptimiser_MCE.pl
#
#        USAGE: ./SSAKEoptimiser_MCE.pl  
#
#  DESCRIPTION: 
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

# code outline
#
# + get input files and sort them by size
# + submit the size sorted files to the SSAKE runs in parallel
# + iterate through the kmers in parallel
# + after kmer iteration for each file,
# 	 - find the output file that contains the longest contig
#	 - remove all other outputfiles
# + print reports of all steps to a log file including processing times
# + first task: creates SSAKE command lines
# + second task: runs SSAKE
# + finds output file with longest contig and remove the rest


use strict;
use warnings;

use MCE::Loop chunk_size => 1;

MCE::Loop::init {
	max_workers => 10
};

my @files = glob("input/*.fa");
#for (@files){ 
#	s/.*\///;
#	print "$_\n"; 
#}

my %h = mce_loop {
	my $file = $_;
	my $wc = `wc -l $_ | cut -f1 -d" "`;
	chomp $wc;
	MCE->gather($file, $wc);
} @files;

#for my $file (sort {$h{$b} <=> $h{$a}} keys %h){  # sorting hash by values, i. e. line counts§
#	my ($stub, $name) = $file =~ m/(.*\/)(.*)/;  # separate the file name from the path
#	print $name, ":		", $h{$stub . $name}, "		lines.\n";
#}

my @files_sorted = sort {$h{$b} <=> $h{$a}} keys %h;
#print "@files_sorted[5 .. $#files_sorted]\n";
#exit;

##########################################################################
## Test MCE::Step
###########################################################################
#use MCE::Step;
#
#sub setup_commandline {
#	my $cmd;
#	MCE->print(MCE->task_wid, ": ", "$_\n");
#	for my $kmer (11..33){
#		$cmd = "SSAKE -kmer " . $kmer . " -in " . $_;
#		MCE->step($cmd, $_);
#	}
#}
#
#sub run_SSAKE {
#	my ($mce, $cmd, $file)  = @_;
#	MCE->print(MCE->task_wid, ": ", $cmd, "\n");
#	MCE->step($cmd, $file);
#}
#
#sub pick_longest_contig {
#	my ($mce, $cmd, $file)  = @_;
#	if($cmd =~ /-kmer 33/){
#		MCE->print(MCE->task_wid, ": ", "found longest contig for file $file\n");
#	}
#}
#
#mce_step {
#	task_name => ['command', 'SSAKE', 'contig'],
#	max_workers => [ 2, 20, 2 ],
#	chunk_size => 1
##}, \&setup_commandline, \&run_SSAKE, \&pick_longest_contig, @files_sorted[5 .. $#files_sorted];
##}, \&setup_commandline, \&run_SSAKE, \&pick_longest_contig, @files_sorted[5 .. 10];
#}, \&setup_commandline, \&run_SSAKE, \&pick_longest_contig, @files_sorted[0..1];
###########################################################################


#use MCE::Flow Sereal => 1;
#use MCE::Queue;
#
#MCE::Flow::init{
#	chunk_size => 1
#}
#
#my $q1 = MCE::Queue->new;
#my $q2 = MCE::Queue->new;
#
#sub task_end {
#	my ($mce, $task_id, $task_name) = @_;
#
#	if( defined $mce->{user_tasks}->[$task_id + 1] ){
#		my $N_workers = $mce->{user_tasks}->[$task_id + 1]->{max_workers};
#
#		if( $task_name eq 'run_SSAKE' ){
#			$q1->enqueue((undef) x ($N_workers) * 2);
#		}elsif( $task_name eq 'contig' ){
#			$q2->enqueue((undef) x ($N_workers) * 2);
#		}
#	}
#
#	return;
#}

##my @input =  @files_sorted[-5 .. -1]; 
##my $input = \@input;
#
#my $input =  [@files_sorted[-5 .. -1]]; 
#foreach my $file (@{$input}){
#	print "$file\n";
#}
#exit;

use Time::HiRes qw(time sleep);
use MCE::Step;

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
	system($cmd) == 0 or die $!;
#	sleep 10;
	MCE->gather($file, $kmer);
}

sub gather_output {
	my %tmp;
	return sub {
		my ($file, $kmer) = @_;
		push @{ $tmp{$file} }, $kmer;
		if( scalar @{ $tmp{$file} } == 33-11+1){
			MCE->print("Finding longest contig for $file.\n");
			find_longest_contig($file);
			delete $tmp{$file};
		}
		return;
	}
	
}

sub find_longest_contig{
	my $file = shift;

	# get output file names containing the contigs
	my @output = `ls $file*contigs`; # use glob
	chomp @output;

	my $longest_contig_length = 0;
	my %longest_contig_hash;

	# get the length of the longest contig for each SSAKE x kmer run
	for my $output (@output){
		if(-z $output){ # if the contig file is empty
			$output =~ s/contigs$/\*/; 
#			MCE->print("Contig file is empty for $output\n");
			system("rm -f $output")==0 or die $!;
			next; 
		}
		my $cmd = "awk \'NR%2==0\' " . $output . " | perl -ne\'chomp; print length, \"\\n\";\' | sort -rn | head -n1";
		$longest_contig_length = `$cmd`;
		chomp $longest_contig_length;
#		MCE->print("$output : $longest_contig_length\n");
		unless ($longest_contig_length){
			MCE->print("Found no contig length for $output\n");
			next;
		}; 
		$longest_contig_hash{$output} = $longest_contig_length;
	}

	# exit if no contigs could be assembled
	unless(keys %longest_contig_hash){
		MCE->print("No contigs could be assembled for $file.\n");
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
	MCE->print("File containing longest contig:\n", "$keep : $longest_contig_hash{$keep}\n\n");


	map { s/contigs$/\*/ } @remove;
	system("rm -f @remove")==0 or die $!;

	return;
}

my $input =  [@files_sorted[30 .. 40]]; 

mce_step {
	task_name => ['command', 'SSAKE'],
	max_workers => [ 1, 20 ],
	chunk_size => 1,
	gather => gather_output()
#}, \&setup_commandline, \&run_SSAKE, \&pick_longest_contig, @files_sorted[5 .. $#files_sorted];
}, \&setup_commandline, \&run_SSAKE, $input;

