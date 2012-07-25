#!/usr/bin/env perl
# purge_PCR_duplicates.pl
use strict; use warnings;
use Getopt::Std;
our ($opt_h, $opt_v, $opt_p, $opt_s);
getopts('hsvp:');
use Data::Dumper;

my $version = "1.0, 22/08/2011, author: Claudius Kerth\n";

my $usage = "
This script takes .fastq input files from standard paired end sRAD libraries (prepared using random shearing). 
It removes PCR duplicates as well as many sequencing errors using the paired end sequence information. It can 
also take the output of the 'process_radtags' quality filtering programme of 'stacks' (Julian Catchen) and 
produce .fastq output files which can then be analysed in a correct way by 'RADtags' of 'RADtools' (John Davey).

usage: ./purge_PCR_duplicates.pl [options] <arguments...>
options:
	-v version
	-h help
	-p <number of processors available>
	-s don't remove PCR duplicates, only output overlapping and matching single end and paired end files
	   (necessary when trying to use 'process_radtags' output with 'RADtags')
	   	
";

die $usage unless @ARGV;

if ($opt_v) { print $version; exit;
}

if ($opt_h) { print $usage; exit;}

my $max_processes = $opt_p;

# Non-core modules
eval {
    require Parallel::ForkManager;
    Parallel::ForkManager->import();
};
die
"RADtags requires the CPAN module Parallel::ForkManager. Please install this package and add it to your Perl library path.\n"
  if $@;
  

#system("date");
my $directory         = '.';

################################################################################
# get a list of filenames to process
################################################################################

# get all .fastq filenames
opendir( my $readdir, $directory )
      or die "Can't open directory $directory!\n";
          my @read_filenames = grep { /(.fastq)$/ } readdir $readdir;
		  my @read_filenames1 = grep { !/^new/ } @read_filenames;
		  @read_filenames = @read_filenames1; 	
    if ( @read_filenames == 0 ) {
        die
"No read files found with extension '.fastq'.\n";
    }
#foreach (@read_filenames) {print; print "\n";} 
closedir $readdir;


my %pair_filenames;
my @paired_files;
foreach my $read_filename (@read_filenames) {
    if ( $read_filename =~ /^(.+)\_1\.fastq$/ ) {
        my $filename_stub = $1;
        foreach my $read2_filename (@read_filenames) {
            # if the single end read file has a corresponding paired end file:
            if ( $read2_filename =~
                /^$filename_stub\_2\.fastq$/ )
            {
                push @paired_files, $read_filename; # store the single end read file name in the array @paired_files
                $pair_filenames{$read_filename} = $read2_filename; # store the paired end read file name in a hash
            }
        }
    }
}
if ( @paired_files > 0 ) { @read_filenames = @paired_files; }

# initialise the parallel processing
my $pm = new Parallel::ForkManager( scalar(@read_filenames) );
$pm->set_max_procs($max_processes);

foreach my $read_filename (@read_filenames) {
	# start this loop for as many filenames as there are processors available
	$pm->start and next;
	
	my ($sample_name) = $read_filename =~ /^(.+)_1\.fastq/;
#   print $sample_name, "\n";

#	open single end file 
    open my $single_in, '<', "$directory/$read_filename"
          or die "Can't open input file $read_filename!\n";
          
#     open paired end file
	my $paired_in;     
    if ( defined $pair_filenames{$read_filename} ) { # if there is a paired end read file for this single end read file     
    	open $paired_in, '<', "$directory/$pair_filenames{$read_filename}"
                  or die "Can't open input file $pair_filenames{$read_filename}!\n";
    } else {
    	next; # get the single end read filename
    }
    
    my $overlap = sort_and_match( $single_in, $paired_in, $sample_name);
    
    	if ($opt_s) 
	{
		print "Created new\_$sample_name\_(1/2)\.fastq files with overlapping and corresponding records.\n";
		$pm->finish;
	}
	
	# store everything in a big hash
	my %uniques = ();
	
	my $paired_end_read_length = weed_out_PCR_dupl(\%uniques, $sample_name);
	
	weed_out_seq_errors(\%uniques, $paired_end_read_length);
	
	print_purged_fastq(\%uniques, $sample_name, $overlap);
	
	$pm->finish;

	
#system("display2fastq.pl -l header $sample_name\_1\.purged $sample_name\_2\.purged | head") == 0 or die "Can't display the two new fastq files: $!";
#	
#system("date");

}
$pm->wait_all_children;

################################################################################
#### Name:       SORT_AND_MATCH
#### Function:   print out sorted and matching single and paired end .fastq files
#### Returns:    nothing (outputs to file)
################################################################################
sub sort_and_match {
	
	my ($single_in, $paired_in, $sample_name) = @_;	
	
#   store the records of the two corresponding single end and paired end files into a big hash
    my %fastq_record;
    
    # read in the first single end record and check format
    my $header    = substr <$single_in>, 0, -2; # remove the last two signs from header (linefeed and '1')
	my $seq       = <$single_in>;
	my $qual_head = <$single_in>;
	my $qual      = <$single_in>;
	die "single end input file not in .fastq format: $sample_name\_1\n" if 
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
		# insert a check for fastq format here
		last if (!$header) || (!$seq) || (!$qual_head) || (!$qual);
		push @{ $fastq_record{paired}{$header} }, ($seq, $qual_head, $qual);
	}
	close $paired_in;
	
	###############################################################################
	# print two new .fastq files sorted by the single end read sequences 
	# and a paired end file that corresponds to the single end file line by line,
	# i. e. print out only single end records which also have a paired end record
	###############################################################################

	open my $new_single, ">", "new\_$sample_name\_1\.fastq" or die "Couldn't open the file new\_$sample_name\_1\.fastq: $!";
	open my $new_paired, ">", "new\_$sample_name\_2\.fastq" or die "Couldn't open the file new\_$sample_name\_2\.fastq: $!";
	
	my $overlap = 0;
	my $single_count = 0;
	foreach my $header ( sort { $fastq_record{single}{$a}[0] cmp $fastq_record{single}{$b}[0] } keys %{ $fastq_record{single} } ){
		if (exists $fastq_record{paired}{$header}) {
			print $new_single 
				$header, "1\n",
				@{ $fastq_record{single}{$header} };
				
			print $new_paired
				$header, "2\n",
				@{ $fastq_record{paired}{$header} };
				
			$overlap++;
		}
		$single_count++;
	}
	print "Found $overlap single end records with corresponding paired end record for sample $sample_name.\n";
	printf "%s%.0f%s", 
		"Discarded ", (1-$overlap/$single_count)*100, "% of single end records because of missing corresponding paired end record (by header).\n";

	close $new_single;
	close $new_paired;
	
	# clear memory
	undef %fastq_record;
	
	return ($overlap);
	
}
################################################################################

################################################################################
##### NAME: weed_out_pcr_dupl
##### FUNCTION: weed out PCR duplicates
##### RETURN: paired end read length
################################################################################

sub weed_out_PCR_dupl {
	
	my ($u, $sample_name) = @_;
	
	# read in the new sorted and overlapping fastq records
	open my $new_single_in, "<", "new\_$sample_name\_1\.fastq" 
		or die "Couldn't open the file new_fastq_single: $!";
	open my $new_paired_in, "<", "new\_$sample_name\_2\.fastq" 
		or die "Couldn't open the file new_fastq_paired: $!";

	
	# initialize $next_single as a global variable so that it can be handed over 
	# from the first while loop structure to the following
	my $next_single; 
	
	
	#######################
	# get the first unique
	#######################
	my $new_unique = get_record( $new_single_in ); 
	# extract the first single end read sequence of the new unique
	my $unique_seq = $new_unique->{seq};
	# push the first single end record of the first unique, 
	# as an anonymous hash on an anonymous array of single end records
	push @{ $u->{$unique_seq}{single} }, $new_unique;
	# push the first paired end record of the first unique, 
	# as an anonymous hash on an anonymous array of paired end records
	push @{ $u->{$unique_seq}{paired} }, get_record ( $new_paired_in );
	$u->{$unique_seq}{paired}[0]{count} = 1;
	my $paired_end_read_length = length($u->{$unique_seq}{paired}[0]{seq}) -1; # don't count the newline character at the end
	$u->{$unique_seq}{PCR_dupl} = 0;
	
	RECORD:
	while ( $next_single = get_record( $new_single_in ) ) { # get the next single end record as an anonymous hash
	#   exit the loop if the next single end read is not equal to the current unique sequence,
	#   i. e. belongs to the next unique
	    if ( $next_single->{seq} ne $unique_seq ) {
	    	my $last_index = $#{ $u->{$unique_seq}{paired} };
			if ( !exists $u->{$unique_seq}{paired}[$last_index]{count} ) { 
	   			$u->{$unique_seq}{paired}[$last_index]{count} = 1;
	    	}
	    	last;
		} 
	#   get the next paired end record, as an anonymous hash
	    my $next_paired = get_record( $new_paired_in );
	    
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
	  				# shortcut the wile loop and get the next single end record
	   	    		next RECORD;
	   		}
	   	}
	#   Add single and paired end records to the 'uniques' hash which only contains
	#   non-PCR duplicate records
	    push @{ $u->{$unique_seq}{single} }, $next_single;
	    push @{ $u->{$unique_seq}{paired} }, $next_paired;
	}
	
	
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
		push @{ $u->{$unique_seq}{paired} }, get_record ( $new_paired_in );
		$u->{$unique_seq}{paired}[0]{count} = 1;
		$u->{$unique_seq}{PCR_dupl} = 0;
		
		RECORD:
		while ( $next_single = get_record( $new_single_in ) ) { # get the next single end record as an anonymous hash
		    if ( $next_single->{seq} ne $unique_seq ) {
		    	my $last_index = $#{ $u->{$unique_seq}{paired} };
				if ( !exists $u->{$unique_seq}{paired}[$last_index]{count} ) { 
	    			$u->{$unique_seq}{paired}[$last_index]{count} = 1;
				}
				last; 
			}
		    my $next_paired = get_record( $new_paired_in );
	    	for my $i ( 0 .. $#{ $u->{$unique_seq}{paired} } ) {
	    		if ( !exists $u->{$unique_seq}{paired}[$i]{count} ) { 
	    			$u->{$unique_seq}{paired}[$i]{count} = 1;
				}
	    		my $dist = ( $next_paired->{seq} ^ $u->{$unique_seq}{paired}[$i]{seq} ) =~ tr/\001-\255//;
		   		if ( $dist < ($paired_end_read_length/5) ){					
	   				$u->{$unique_seq}{paired}[$i]{count}++; 
	   				$u->{$unique_seq}{PCR_dupl}++;
	    	    	next RECORD;
		   		}
	    	}
		    push @{ $u->{$unique_seq}{single} }, $next_single;
		    push @{ $u->{$unique_seq}{paired} }, $next_paired;
		}
	
	}
	
	system("rm new_$sample_name\_?.fastq") == 0 or die "Can't remove intermediate files: $!"; 
	
	return ($paired_end_read_length);
}
################################################################################


################################################################################
##### NAME: weed_out_seq_errors
##### FUNCTION: iterate over the new hash and weed out sequencing errors
##### RETURNS: nothing, modifies the global %uniques
################################################################################

sub weed_out_seq_errors {
	
	my ($u, $paired_end_read_length) = @_;
		
#	 get a sorted list if uniques
	my @uniques_sorted = sort keys %{ $u };
	
#	iterate over uniques, deleting them if they are sequencing errors:
	my $i = 0;
	while ( $i <= $#uniques_sorted -1) { # the last index of the list of sorted uniques changes when a unique is deleted
		my $uniq_seq = $uniques_sorted[$i]; # get one unique sequence from the sorted list
#		print $uniq_seq;
		my $next_uniq = $uniques_sorted[$i+1]; # get the next unique sequence from the sorted list
#		print $next_uniq;
#		 iterate over the different paired end sequences of the first unique, 
#		 deleting it and its single end record if its single end read contains 
#		 a sequencing error:
		my $j = 0;
	PAIRED1:	
		while  ( $j <= $#{ $u->{$uniq_seq}{paired} } ){
#			 iterate over paired end sequences of the second unique,
#			 deleting them together with their single end records if their single end read contain an error
			my $k = 0;
			while ( $k <= $#{ $u->{$next_uniq}{paired} } ){
				my $dist = 
				( $u->{$uniq_seq}{paired}[$j]{seq} ^ $u->{$next_uniq}{paired}[$k]{seq} )
				=~ tr/\001-\255//;
#				 if the two paired end sequences are the same
				if ( $dist < ($paired_end_read_length/5) ){
#					 if the paired end read from the first unique has only a count of 1 but the other a higher count
					if ( $u->{$uniq_seq}{paired}[$j]{count} == 1 and $u->{$next_uniq}{paired}[$k]{count} > 1 ){
#						 remove single and paired end records of the erroneous read from the hash
#						print "removed record number ", $j+1, " for first unique: $uniques{$uniq_seq}{paired}[$j]{seq}";
						splice @{ $u->{$uniq_seq}{paired} }, $j;
						splice @{ $u->{$uniq_seq}{single} }, $j;
#						 go to the next paired end record of the first unique without incrementing $j
#						 because the array @{ $uniques{$uniq_seq}{paired} } has just lost one element:
						next PAIRED1;
					}
#					 if the paired end read of the next unique has only a count of 1 but the other a higher count	 
					elsif ( $u->{$uniq_seq}{paired}[$j]{count} > 1 and $u->{$next_uniq}{paired}[$k]{count} == 1 ){
#						 remove single and paired end records of the erroneous read from the hash
#						print "removed record number ", $k+1, " for next unique: $uniques{$next_uniq}{paired}[$k]{seq}";
						splice @{ $u->{$next_uniq}{paired} }, $k;
						splice @{ $u->{$next_uniq}{single} }, $k;
#						 go to the next paired end record of the next unique without incrementing $k
#						 because the array @{ $uniques{$next_uniq}{paired} } has just lost an element:
						next;
					}
#					 if both paired end reads have a count of 1, I cannot decide which record is the one with the sequencing error:
					else { 
						$k++; # increment the index for the paired end records of the 'next_uniq'
						next;
					}
				}
#				 if both paired end sequences are different, increment the index of the paired-end records of the 'next_uniq'
				else { $k++;}	
			}
#			 once you have been through all paired end records of the 'next_uniq', 
#			 increment the index for the paired end records of the 'uniq_seq':
			$j++; 
		}
#		 if the first unique sequence and the next unique sequence now still have records associated with them:
		if ( @{ $u->{$uniq_seq}{single} } and @{ $u->{$next_uniq}{single} } ) { 
			$i++; # increment the index for the next unique sequence in @uniques_sorted
		}
#		 if the first unique sequence now has no more records associated with it 
#		 (because it contained only sequencing errors whose records got deleted):
		elsif ( !@{ $u->{$uniq_seq}{single} } ){ 
#			 delete the unique
			delete $u->{$uniq_seq};
			splice @uniques_sorted, $i, 1;
		}
#		 if the next unique sequence now has no more records associated with it
#		 (because it contained only sequencing errors whose records got deleted):
		elsif ( !@{ $u->{$next_uniq}{single} } ){ 
#			 delete the unique
			delete $u->{$next_uniq};
			splice @uniques_sorted, $i+1, 1;
		}
		else { print "Error!\n";}
	}
	
	return;
	
}
################################################################################


###############################################################################
##### NAME: print_purged_fastq
##### FUNCTION: print out result
##### RETURN: nathing (prints to output)
#############################################################################
	
sub print_purged_fastq {
	
	my ($u, $sample_name, $overlap) = @_;
		
	
	open my $new_single_out, ">", "$sample_name\_1.purged" or die "Can't open file $sample_name\_1.purged: $!";
	open my $new_paired_out, ">", "$sample_name\_2.purged" or die "Can't open file $sample_name\_2.purged: $!";
	
	print "Writing to $sample_name\_1.purged and $sample_name\_2.purged output files ... \n";

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
	printf "%s%.0f%s", "Removed ", (1-$count/$overlap)*100, "% of all overlapping records, which were either PCR duplicates or contained sequencing errors.\n";

	return;
	
}



################################################################################
# NAME: get_record
# GETS: filehandle
# RETURNS: one record in an anonymous hash
################################################################################

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
################################################################################
