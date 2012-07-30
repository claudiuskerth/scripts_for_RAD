#!/usr/bin/env perl
# purge_PCR_duplicates.pl
use strict; use warnings;

use Data::Dumper;

my $version = "1.1, 30/07/2012, author: Claudius Kerth\n";

 # :TODO      :07/29/2012 17:54:53:CEK: avoid writing an intermediate file to disk with the subroutine sort_and_match
 # :TODO      :07/29/2012 17:56:08:CEK: keep a consensus sequence of paired end sequences from a set of PCR duplicates
 # :TODO      :07/29/2012 17:57:14:CEK: correct obvious sequencing errors instead of throwing error containing reads away

my $usage = "
This script takes .fastq input files from standard paired end sRAD libraries (prepared using random shearing) from the
current directory. It removes PCR duplicates while keeping the median of quality scores for the single-end (SE) reads 
that belong to a fragment (i.e. all but one of which are PCR duplicates). PCR duplicates are identified by almost identical 
paired-end (PE) read sequences (i. e. allowing for sequencing errors). The script currently writes an intermediate file 
(for each input) to disk in the current directory, which should be the same size as in the input file. Input files must 
follow a strict naming convention:
a set of paired files (i. e. SE and corresponding PE file) must have the same file name stub and end with either
\"fq_1\" for the SE file or \"fq_2\" for the PE file. Example: XXX.fq_1 and XXX.fq_2 . Also, the script assumes that the headers
of the fastq records have not been manipulated. The script is parallelized to speed up processing. It is highly recommended to 
run this script on demultiplexed sequence files and not on the whole data of a HiSeq lane, for instance, unless you have a lot
of memory in your computer. The output files contain the keyword \"purged\" in their names.

usage: ./purge_PCR_duplicates.pl [options] 

options:
	-v version
	-h help
	-p <number of processors available> (default: 1)
	   	
";

my $max_processes = 1;

parse_command_line();

#===  FUNCTION  ================================================================
#         NAME: PARSE_COMMAND_LINE
#      PURPOSE: parses @ARGV for parameter input on the command line
#   PARAMETERS: ????
#      RETURNS: ????
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub parse_command_line{
#	die $usage unless @ARGV;
	while(@ARGV){
		$_ = shift @ARGV;
		if (/^-v$/) {die $version;}
		if (/^-h|-help$/) {die $usage;}
		if (/^-p$/) {$max_processes = shift @ARGV;}
	}
} ## --- end sub PARSE_COMMAND_LINE


# Non-core modules
eval {
    require Parallel::ForkManager;
    Parallel::ForkManager->import();
};
die
"purge_PCR_duplicates.pl requires the CPAN module Parallel::ForkManager. Please install this package and add it to your Perl library path.\n"
  if $@;
  

#system("date");
my $directory         = '.';

################################################################################
# get a list of filenames to process
################################################################################

# get all .fastq filenames
opendir( my $readdir, $directory ) or die "Can't open directory $directory!\n";
my @read_filenames = grep { /(.fq_[1|2])$/ && !/_sorted\.fq/ && !/_purged\.fq/ } readdir $readdir;
if ( @read_filenames == 0 ) { die "\nNo read files found with extension '.fq_1'.\n\n"; }
#foreach (@read_filenames) {print; print "\n";} 
closedir $readdir;
#exit;

my %pair_filenames;
my @paired_files; # shall contain the SE filenames, which have a corresponding PE filename
my $filename_stub;

foreach my $SE_filename (grep {/fq_1$/} @read_filenames) {
	($filename_stub) = $SE_filename =~ /^(.+)\.fq_1$/; 
#	print $filename_stub, "\n";
	foreach my $PE_filename (grep {/fq_2$/} @read_filenames) {
#		print $PE_filename, "\n";
		# if the single end read file has a corresponding paired end file:
		if ( $PE_filename =~ /^$filename_stub\.fq_2$/ ) {
			push @paired_files, $SE_filename; # store the single end read file name in the array @paired_files
			$pair_filenames{$SE_filename} = $PE_filename; # store the paired end read file name in a hash
		}
	}
}
#foreach my $SE (sort keys %pair_filenames){
#	print $SE, "\t", "$pair_filenames{$SE}\n";
#}

unless (@paired_files > 0) { die "\nNo paired-end sequence files found! 
Make sure that names of paired files are identical except for their file endings \"fq_1\" or \"fq_2\".\n\n"; }
#

# initialise the parallel processing
my $pm = new Parallel::ForkManager( scalar(@read_filenames) );
$pm->set_max_procs($max_processes);

foreach my $filename (@paired_files) {
	# start this loop for as many filenames as there are processors available
	$pm->start and next;
	
	my ($sample_name) = $filename =~ /^(.+)\.fq_1$/;

#	open single end file 
    open my $single_in, '<', "$directory/$filename"
		or die "Can't open input file $filename!\n";
          
#   open paired end file
	my $paired_in;     
   	open $paired_in, '<', "$directory/$pair_filenames{$filename}"
		or die "Can't open input file $pair_filenames{$filename}!\n";
    
	# earlier versions of stacks' "process_radtags" filtered single-end reads 
	# without deleting their corresponding paired-end reads. The following subroutine
	# tries to fix that. However with newer versions of "process_radtags" that is
	# not necessary anymore.
	my $overlap = 0;
    $overlap = sort_and_match( $single_in, $paired_in, $sample_name );
    
#	# store everything in a big hash
	my %uniques = ();
	
	weed_out_PCR_dupl(\%uniques, $sample_name, $single_in, $paired_in);
	
	# under further construction
##	weed_out_seq_errors(\%uniques, $paired_end_read_length);
	
	print_purged_fastq(\%uniques, $sample_name, $overlap);
	
	$pm->finish;

	
}
$pm->wait_all_children;
exit;


################################################################################
##### Name:       SORT_AND_MATCH
##### Function:   print out sorted and matching single and paired end fastq files
##### Returns:    nothing (outputs to file)
#################################################################################
sub sort_and_match {
	
	my ($single_in, $paired_in, $sample_name) = @_;	
	
#   store the records of the two corresponding single end and paired end files into a big hash
    my %fastq_record;
    
    # read in the first single end record and check format
    my $header    = substr <$single_in>, 0, -2; # remove the last two signs from header (linefeed and '1')
	my $seq       = <$single_in>;
	my $qual_head = <$single_in>;
	my $qual      = <$single_in>;
	die "single end input file not in .fastq format: $sample_name\n" if 
		(
		$header !~ /^@.+$/ 
		or $seq =~ /^[^AGCTN]$/
		or $qual_head !~ /^\+/
		or $qual =~ /^[^!-h]$/
		);
	last if (!$header) || (!$seq) || (!$qual_head) || (!$qual);
	push @{ $fastq_record{single}{$header} }, ($seq, $qual_head, $qual);

    # read in the rest of the single end records
    while (<$single_in>) {
		$header    = substr $_, 0, -2; # remove the last two signs from header (linefeed and '1')
		$seq       = <$single_in>;
		$qual_head = <$single_in>;
		$qual      = <$single_in>;
		last if (!$header) || (!$seq) || (!$qual_head) || (!$qual);
		push @{ $fastq_record{single}{$header} }, ($seq, $qual_head, $qual);
	}
	close $single_in;
	
	# read in the first paired end record and check format
	$header    = substr <$paired_in>, 0, -2;
	$seq       = <$paired_in>;
	$qual_head = <$paired_in>;
	$qual      = <$paired_in>;
	die "paired end input file not in .fastq format: $sample_name\_2\n" if 
		(
		$header !~ /^@.+$/ 
		or $seq =~ /^[^AGCTN]$/
		or $qual_head !~ /^\+/
		or $qual =~ /^[^!-h]$/
		);
	last if (!$header) || (!$seq) || (!$qual_head) || (!$qual);		
	push @{ $fastq_record{paired}{$header} }, ($seq, $qual_head, $qual);

	# read in the rest of the paired end records
	while (<$paired_in>) {
		$header    = substr $_, 0, -2;
		$seq       = <$paired_in>;
		$qual_head = <$paired_in>;
		$qual      = <$paired_in>;
		last if (!$header) || (!$seq) || (!$qual_head) || (!$qual);
		push @{ $fastq_record{paired}{$header} }, ($seq, $qual_head, $qual);
	}
	close $paired_in;
	
	###############################################################################
	# print two new fastq files sorted by the single end read sequences 
	# and a paired end file that corresponds to the single end file line by line,
	# i. e. print out only single end records which also have a paired end record
	###############################################################################

	open (my $SE_sorted, ">", "$sample_name\_sorted.fq_1") or die "Couldn't open the file $sample_name\_sorted.fq_1: $!";
	open (my $PE_sorted, ">", "$sample_name\_sorted.fq_2") or die "Couldn't open the file $sample_name\_sorted.fq_2: $!";
	
	my $overlap = 0;
	my $single_count = 0;
	

	# sort SE records lexically on sequence
	foreach my $header ( sort { $fastq_record{single}{$a}[0] cmp $fastq_record{single}{$b}[0] } keys %{ $fastq_record{single} } ){
		# the corresponding paired-end record is found by its header
		if (exists $fastq_record{paired}{$header}) {
			print $SE_sorted
				$header, "1\n",
				@{ $fastq_record{single}{$header} };
				
			print $PE_sorted
				$header, "2\n",
				@{ $fastq_record{paired}{$header} };
				
			$overlap++;
		}
		$single_count++;
	}
	print "Found $overlap single end records with corresponding paired end record for sample $sample_name.\n";
	printf "%s%.0f%s", 
		"Discarded ", (1-$overlap/$single_count)*100, "% of single end records because of missing corresponding paired end record (by header).\n";

	close $SE_sorted;
	close $PE_sorted;
	
	# clear memory
	undef %fastq_record;

	return($overlap);
	
} # -- end of sort_and_match
#################################################################################

#################################################################################
###### NAME: weed_out_pcr_dupl
###### FUNCTION: weed out PCR duplicates
###### RETURN: paired end read length
#################################################################################

sub weed_out_PCR_dupl {

	my ($u, $sample_name) = @_;

	# read in the new sorted and overlapping fastq records
	open my $sorted_single_in, "<", "$sample_name\_sorted.fq_1" 
		or die "Couldn't open the file $sample_name\_sorted.fq_1: $!";
	open my $sorted_paired_in, "<", "$sample_name\_sorted.fq_2" 
		or die "Couldn't open the $sample_name\_sorted.fq_2: $!";

	# initialize $next_single as a global variable so that it can be handed over 
	# from the first while loop structure to the following
	my $next_single; 

	#######################
	# get the first unique
	#######################
	my $new_unique = get_record( $sorted_single_in ); 
	# extract the first single end read sequence of the new unique
	my $unique_seq = $new_unique->{seq};
	# push the first single end record of the first unique, 
	# as an anonymous hash on an anonymous array of single end records:
	push @{ $u->{$unique_seq}{single} }, $new_unique;
	# push the first paired end record of the first unique, 
	# as an anonymous hash on an anonymous array of paired end records
	push @{ $u->{$unique_seq}{paired} }, get_record ( $sorted_paired_in );
	$u->{$unique_seq}{paired}[0]{count} = 1;
	my $paired_end_read_length = length($u->{$unique_seq}{paired}[0]{seq}) -1; # don't count the newline character at the end
	$u->{$unique_seq}{PCR_dupl} = 0;
	# I want to store the median of each quality score position for a set of PCR duplicates. A set of PCR
	# duplicates is defined by identical SE and PE sequence.
	push @{ $u->{$unique_seq}{paired}[0]{SE_qual} }, $new_unique->{qual};


	RECORD:
	while ( $next_single = get_record( $sorted_single_in ) ) { # get the next single end record as an anonymous hash
	#   exit the loop if the next single end read is not equal to the current unique sequence,
	#   i. e. belongs to the next unique
		if ( $next_single->{seq} ne $unique_seq ) {
			my $last_index = $#{ $u->{$unique_seq}{paired} };
			if ( !exists $u->{$unique_seq}{paired}[$last_index]{count} ) { 
				$u->{$unique_seq}{paired}[$last_index]{count} = 1;
			}

			calculate_median_qual($u, $unique_seq);

			last;
		} 
	#   get the next paired end record, as an anonymous hash
		my $next_paired = get_record( $sorted_paired_in );
		
	#   Check if the paired end sequence has already been seen for this unique.
	#	If so, skip recording both single and paired end records because they belong to a PCR duplicate.

	#	Iterate over all previously recorded paired end sequences for that unique.
		for my $i ( 0 .. $#{ $u->{$unique_seq}{paired} } ) {
			# initalize the count for each paired end sequence to 1, if it not already exists 
			if ( !exists $u->{$unique_seq}{paired}[$i]{count} ) { 
				$u->{$unique_seq}{paired}[$i]{count} = 1;
			}
			# count bp mismatches, see http://www.perlmonks.org/?node_id=500235 for explanation
			my $dist = ( $next_paired->{seq} ^ $u->{$unique_seq}{paired}[$i]{seq} ) =~ tr/\001-\255//;

			# If the distance between the two paired end sequences is less than 1/5 of their length ... .
			# By chance, two random sequences should have on average 1/4 of their bp matching.
			# If paired end read length is 51, then a dist of 10 will still be regarded as equal, 
			# i. e. 10 PCR or sequencing errors in either paired end sequence.				
			if ( $dist < ($paired_end_read_length/5) ) {
				# increment read count for each fragment (not fragment count as in RADtags!)
				$u->{$unique_seq}{paired}[$i]{count}++; 
				# increment total count of PCR duplicates for the unique
				$u->{$unique_seq}{PCR_dupl}++; 
				# keep the SE quality scores for this PCR duplicate:
				push @{ $u->{$unique_seq}{paired}[$i]{SE_qual} }, $next_single->{qual};
				# shortcut the while loop and get the next single end record
				next RECORD;
			}
		}
	#   Add single and paired end records to the 'uniques' hash which only contains
	#   non-PCR duplicate records
		push @{ $u->{$unique_seq}{single} }, $next_single;
		push @{ $u->{$unique_seq}{paired} }, $next_paired;
		# get the currently last index of the array of paired-end records for this SE unique
		my $last_index = $#{ $u->{$unique_seq}{paired} };
		# store the SE quality scores of this first read from a possible set of PCR duplicates 
		push @{ $u->{$unique_seq}{paired}[$last_index]{SE_qual} }, $next_single->{qual}; 
	}

#	print Dumper($u->{$unique_seq}{paired});
#	print "\n***\n";
#	print Dumper($u->{$unique_seq}{single});
#	
#	return;


	################################
	# now get all the other uniques
	################################

	while ( my $new_unique = $next_single ) {
		# extract the first single end read sequence of the new unique
		my $unique_seq = $new_unique->{seq};
		# push the first single end record of the first unique, 
		# as an anonymous hash on an anonymous array of single end records
		push @{ $u->{$unique_seq}{single} }, $new_unique;
		# push the first paired end record of the first unique, 
		# as an anonymous hash on an anonymous array of paired end records
		push @{ $u->{$unique_seq}{paired} }, get_record ( $sorted_paired_in );
		$u->{$unique_seq}{paired}[0]{count} = 1;
		$u->{$unique_seq}{PCR_dupl} = 0;
		push @{ $u->{$unique_seq}{paired}[0]{SE_qual} }, $new_unique->{qual};
		
		RECORD:
		while ( $next_single = get_record( $sorted_single_in ) ) { # get the next single end record as an anonymous hash
			if ( $next_single->{seq} ne $unique_seq ) {
				my $last_index = $#{ $u->{$unique_seq}{paired} };
				if ( !exists $u->{$unique_seq}{paired}[$last_index]{count} ) { 
					$u->{$unique_seq}{paired}[$last_index]{count} = 1;
				}

				calculate_median_qual($u, $unique_seq);

				last; 
			}
			my $next_paired = get_record( $sorted_paired_in );
			for my $i ( 0 .. $#{ $u->{$unique_seq}{paired} } ) {
				if ( !exists $u->{$unique_seq}{paired}[$i]{count} ) { 
					$u->{$unique_seq}{paired}[$i]{count} = 1;
				}
				my $dist = ( $next_paired->{seq} ^ $u->{$unique_seq}{paired}[$i]{seq} ) =~ tr/\001-\255//;
				if ( $dist < ($paired_end_read_length/5) ){					
					$u->{$unique_seq}{paired}[$i]{count}++; 
					$u->{$unique_seq}{PCR_dupl}++;
					# keep the SE quality scores for this PCR duplicate:
					push @{ $u->{$unique_seq}{paired}[$i]{SE_qual} }, $next_single->{qual};
					next RECORD;
				}
			}
			# if the new record is not a PCR duplicate of a record already seen
			# then add it to the hash of unique fragments
			push @{ $u->{$unique_seq}{single} }, $next_single;
			push @{ $u->{$unique_seq}{paired} }, $next_paired;
			# get the currently last index of the array of paired-end records for this SE unique
			my $last_index = $#{ $u->{$unique_seq}{paired} };
			# store the SE quality scores of this first read from a possible set of PCR duplicates 
			push @{ $u->{$unique_seq}{paired}[$last_index]{SE_qual} }, $next_single->{qual}; 
		}

	}

	system("rm $sample_name\_sorted.fq_?") == 0 or die "Can't remove intermediate files: $!"; 

#	print Dumper($u);

} # -- end of weed_out_PCR_dupl

#################################################################################
#
#
#################################################################################
###### NAME: weed_out_seq_errors
###### FUNCTION: iterate over the new hash and weed out sequencing errors
###### RETURNS: nothing, modifies the global %uniques
#################################################################################
#
#sub weed_out_seq_errors {
#	
#	my ($u, $paired_end_read_length) = @_;
#		
##	 get a sorted list if uniques
#	my @uniques_sorted = sort keys %{ $u };
#	
##	iterate over uniques, deleting them if they are sequencing errors:
#	my $i = 0;
#	while ( $i <= $#uniques_sorted -1) { # the last index of the list of sorted uniques changes when a unique is deleted
#		my $uniq_seq = $uniques_sorted[$i]; # get one unique sequence from the sorted list
##		print $uniq_seq;
#		my $next_uniq = $uniques_sorted[$i+1]; # get the next unique sequence from the sorted list
##		print $next_uniq;
##		 iterate over the different paired end sequences of the first unique, 
##		 deleting it and its single end record if its single end read contains 
##		 a sequencing error:
#		my $j = 0;
#	PAIRED1:	
#		while  ( $j <= $#{ $u->{$uniq_seq}{paired} } ){
##			 iterate over paired end sequences of the second unique,
##			 deleting them together with their single end records if their single end read contain an error
#			my $k = 0;
#			while ( $k <= $#{ $u->{$next_uniq}{paired} } ){
#				my $dist = 
#				( $u->{$uniq_seq}{paired}[$j]{seq} ^ $u->{$next_uniq}{paired}[$k]{seq} )
#				=~ tr/\001-\255//;
##				 if the two paired end sequences are the same
#				if ( $dist < ($paired_end_read_length/5) ){
##					 if the paired end read from the first unique has only a count of 1 but the other a higher count
#					if ( $u->{$uniq_seq}{paired}[$j]{count} == 1 and $u->{$next_uniq}{paired}[$k]{count} > 1 ){
##						 remove single and paired end records of the erroneous read from the hash
##						print "removed record number ", $j+1, " for first unique: $uniques{$uniq_seq}{paired}[$j]{seq}";
#						splice @{ $u->{$uniq_seq}{paired} }, $j;
#						splice @{ $u->{$uniq_seq}{single} }, $j;
##						 go to the next paired end record of the first unique without incrementing $j
##						 because the array @{ $uniques{$uniq_seq}{paired} } has just lost one element:
#						next PAIRED1;
#					}
##					 if the paired end read of the next unique has only a count of 1 but the other a higher count	 
#					elsif ( $u->{$uniq_seq}{paired}[$j]{count} > 1 and $u->{$next_uniq}{paired}[$k]{count} == 1 ){
##						 remove single and paired end records of the erroneous read from the hash
##						print "removed record number ", $k+1, " for next unique: $uniques{$next_uniq}{paired}[$k]{seq}";
#						splice @{ $u->{$next_uniq}{paired} }, $k;
#						splice @{ $u->{$next_uniq}{single} }, $k;
##						 go to the next paired end record of the next unique without incrementing $k
##						 because the array @{ $uniques{$next_uniq}{paired} } has just lost an element:
#						next;
#					}
##					 if both paired end reads have a count of 1, I cannot decide which record is the one with the sequencing error:
#					else { 
#						$k++; # increment the index for the paired end records of the 'next_uniq'
#						next;
#					}
#				}
##				 if both paired end sequences are different, increment the index of the paired-end records of the 'next_uniq'
#				else { $k++;}	
#			}
##			 once you have been through all paired end records of the 'next_uniq', 
##			 increment the index for the paired end records of the 'uniq_seq':
#			$j++; 
#		}
##		 if the first unique sequence and the next unique sequence now still have records associated with them:
#		if ( @{ $u->{$uniq_seq}{single} } and @{ $u->{$next_uniq}{single} } ) { 
#			$i++; # increment the index for the next unique sequence in @uniques_sorted
#		}
##		 if the first unique sequence now has no more records associated with it 
##		 (because it contained only sequencing errors whose records got deleted):
#		elsif ( !@{ $u->{$uniq_seq}{single} } ){ 
##			 delete the unique
#			delete $u->{$uniq_seq};
#			splice @uniques_sorted, $i, 1;
#		}
##		 if the next unique sequence now has no more records associated with it
##		 (because it contained only sequencing errors whose records got deleted):
#		elsif ( !@{ $u->{$next_uniq}{single} } ){ 
##			 delete the unique
#			delete $u->{$next_uniq};
#			splice @uniques_sorted, $i+1, 1;
#		}
#		else { print "Error!\n";}
#	}
#	
#	return;
#	
#}
#################################################################################
#
#
################################################################################
###### NAME: print_purged_fastq
###### FUNCTION: print out result
###### RETURN: nothing (prints to output)
##############################################################################
	
sub print_purged_fastq {
	
#	my ($u, $sample_name, $overlap) = @_;
	my $u = shift @_;
	my $sample_name = shift @_;
	my $overlap = shift @_;
		
	open my $new_single_out, ">", "$sample_name\_purged.fq_1" or die "Can't open file $sample_name\_purged.fq_1: $!";
	open my $new_paired_out, ">", "$sample_name\_purged.fq_2" or die "Can't open file $sample_name\_purged.fq_2: $!";
	
	print "Writing to $sample_name\_purged.fq_2 and $sample_name\_purged.fq_2 output files ... \n";

	my $count = 0;
	foreach my $unique_seq ( sort keys %{ $u }) {
#	print $unique_seq, "\n";	
		for my $i (0 .. $#{ $u->{$unique_seq}{single} } ){
#			chomp $u->{$unique_seq}{single}[$i]{seq};
			print $new_single_out
			$u->{$unique_seq}{single}[$i]{header},
			$u->{$unique_seq}{single}[$i]{seq},
#			" ", $u->{$unique_seq}{PCR_dupl}, "\n",
			$u->{$unique_seq}{single}[$i]{qual_head},
			$u->{$unique_seq}{single}[$i]{qual};
#			
#			chomp $u->{$unique_seq}{paired}[$i]{seq};
			print $new_paired_out
			$u->{$unique_seq}{paired}[$i]{header},
			$u->{$unique_seq}{paired}[$i]{seq},
#			" ", $u->{$unique_seq}{paired}[$i]{count}, "\n",
			$u->{$unique_seq}{paired}[$i]{qual_head},
			$u->{$unique_seq}{paired}[$i]{qual};
			
			$count++;		
		}
	}
	print "The purged files for $sample_name contain each $count records.\n";
	printf "%s%.0f%s", "Removed ", (1-$count/$overlap)*100, 
			"% of all overlapping records, which were either PCR duplicates or contained sequencing errors.\n";

	return;
	
}# -- end of print_purged_fastq



#################################################################################
## NAME: get_record
## GETS: filehandle
## RETURNS: one record in an anonymous hash
#################################################################################

sub get_record {
	my ($fh)  = @_;

	my $header    = <$fh>;
	my $seq       = <$fh>;
	my $qual_head = <$fh>;
	my $qual      = <$fh>;
	
	if ( (!$header) || (!$seq) || (!$qual_head) || (!$qual) ) {
		return;}
	else{
#		print $header, $seq, $qual_head, $qual;
		# return the record as an anonymous hash
		return 
		{
			header => $header,
			seq => $seq,
			qual_head => $qual_head,
			qual => $qual,
		};
	}

}
#################################################################################


#===  FUNCTION  ================================================================
#         NAME: calculate_median_qual
#      PURPOSE: calculates the median of quality scores for all SE reads that
#      			belong to a unique fragment
#   PARAMETERS: reference to the uniques hash, the current unique seq
#      RETURNS: nothing, modifies uniques hash
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub calculate_median_qual {
	my ($u, $unique_seq) = @_;

	# foreach fragment (i. e. unique in SE and PE sequence)
	for (my $i = 0; $i < @{ $u->{$unique_seq}{paired} }; $i++){
		# get the array of quality scores for each SE read that belongs to that fragment
		my @qual_array = @{ $u->{$unique_seq}{paired}[$i]{SE_qual} };
		# foreach position in the quality string
		for (my $j = 0; $j < length($unique_seq); $j++){
			my @median_qual = ();
			# foreach SE read that belongs to that fragment
			foreach my $qual (@qual_array){
				# store the quality scores for that position in an array
				push @median_qual, ord( substr( $qual, $j, 1) );
			}
			# sort quality scores ascending
			@median_qual = sort {$a <=> $b} @median_qual;
			my $index = @median_qual/2; 
			my $median = $median_qual[$index]; # for an uneven length of @median_qual this
			# gives the exact median; for an even length of @median_qual this will give
			# the higher of the two quality scores around the (x.5) median. 

			# now, replace the quality score (earlier taken from the first SE read of that fragment)
			# with the median quality score for all SE reads of that fragment (for the
			# current position in the quality string) 
			substr( $u->{$unique_seq}{single}[$i]{qual}, $j, 1, chr($median) );
		}
	}	
}
