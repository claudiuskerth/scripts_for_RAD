#!/usr/bin/env perl
use strict; use warnings;

#system("rm *trueTags*");

# read in paired file names
# for each pair of .fastq files: read in records, search sequences for remainder of restriction 
# site sequence at the beginning, print out record if true RAD tag sequence

my $usage = "
\nThis script takes a number of paired .fastq file names and removes records if one of their reads
does not begin with the expected remainder of the restriction site sequence.
The input files must be in fastq format. Quality scores in Phred+64 and Phred+33 are accepted.
Paired files must have corresponding records in the same order and file names must be the same,
except for the file name extension (which should be preceded by the only dot in the file name).

$0

-s   single end file name or wildcard expression like \"*.fq_1\"
     (wildcards must be protected with quotation marks against shell expansion)
-p   paired end file name or wildcard expression like \"*.fq_2\"
-f1  pattern to search for in the single end sequences 
     Can include perl regular expressions (e. g. \"^CGAG\", which searches for \"CGAG\" 
     at the beginning of the sequences) when search for exact matches is switched on. 
     Otherwise, the given pattern must only contain the string to search for.
-f2  pattern to search for in the paired end sequences. Comments for \"f1\" apply.
-t   number of processors available
-e   filter only exact matches. By default, the script searches for the given pattern at the 
     beginning of sequences and allows for one mismatch
-m   number of mismatches allowed when filtering for sequences with matching patterns
     Is ignored when exact matches is switched on.
-pos starting position of the search sub-sequence in the read, starting with 1 (default is 1)
     has only an effect on non-exact searches, i. e. when \"-e\" switch is not turned on
-inv invert match, print out records that don't match the search pattern
\n";

#print "\n****printing to output files is currently outcommented for testing purposes!****\n";

print STDERR `date`;

# Non-core modules
eval {
     require Parallel::ForkManager;
     Parallel::ForkManager->import();
     };
die "grep_true_RADtags.pl requires the CPAN module Parallel::ForkManager for parallelising. Please install this package and add it to your Perl library path.\n" if $@;

# initialize variables
my $max_processes = 1;
my $mismatch = 1;
my $start_position = 0;
my ($glob_single, $glob_paired, $pattern_s, $pattern_p, $filtered, $exact, $invert_match);
my @arg = @ARGV;


parse_command_line();

unless ( defined($glob_single) && defined($glob_paired) && ( defined($pattern_s) or defined($pattern_p) ) ) {
	print "\nyour command line:\n", $0, " ",  "@arg\n";
        print $usage;
        exit;
}

my @single_end_files = sort glob($glob_single);
my @paired_end_files = sort glob($glob_paired);

# if no single-end or paired-end files are found 
unless ( @single_end_files and @paired_end_files ) { 
	print "\nyour command line:\n", $0, " ",  "@arg\n";
	print $usage; 
	exit; 
}

print STDERR "command line used for this run:\n", $0, " ", "@arg", $usage;

# initialize parallel processing
my $pm = new Parallel::ForkManager( scalar(@single_end_files) );
$pm->set_max_procs($max_processes);

my ($discarded, $kept, $fname_s_stub);

# iterate over the paired files
foreach my $file ( 0..$#single_end_files ) {
       print "Scanning reads from files: ", $single_end_files[$file], "\t", $paired_end_files[$file], "\n";

       # start this loop for as many filenames as there are processors available
       $pm->start and next;

       open (my $single_in, "<", $single_end_files[$file]) or die "Can't open " . $single_end_files[$file] . "$!\n";
       open (my $paired_in, "<", $paired_end_files[$file]) or die "Can't open " . $paired_end_files[$file] . "$!\n";

       ($discarded, $kept, $fname_s_stub) = filter_true_tags (
		{
			fh_single => $single_in,
			fh_paired => $paired_in,
			single_file_name => $single_end_files[$file],
			paired_end_name => $paired_end_files[$file]
		}
       );	 
       printf "%s%d%s%.2f%s", 
       		"discarded ", $discarded, " records from sample \"$fname_s_stub\" (", $discarded*100/$kept, "%)\n"; 
       printf STDERR "%s%d%s%.2f%s", 
       		"discarded ", $discarded, " records from sample \"$fname_s_stub\" (", $discarded*100/$kept, "%)\n"; 
       $pm->finish;       
}
$pm->wait_all_children;
print "Finished\n";
print STDERR `date`;
exit;

#############################################################################
### FUNCTION: searches sequences in paired files for specified patterns
###	      and prints only those records where the pattern has been found 
### RECEIVES: file handles and file names of paired files to process
### RETURNS:  percentage discarded records, sample name
#############################################################################
sub filter_true_tags {

	my $arg_ref = shift;

	# initialize variables
	my $fh_single = $arg_ref->{fh_single};
	my $fh_paired = $arg_ref->{fh_paired};
	my $fname_s = $arg_ref->{single_file_name};
	my $fname_p = $arg_ref->{paired_end_name};
	my ($fname_s_stub, $fname_p_stub);
	my $discarded = 0;
	my $kept = 0;

	# get sample names
	unless ( ($fname_s_stub ) = ($fname_s =~ /(.+)\..+/) ) { die $usage, $!; }	
	unless ( ($fname_p_stub ) = ($fname_p =~ /(.+)\..+/) ) { die $usage, $!; }
	
	if ( $fname_s_stub ne $fname_p_stub ) { die "paired files don't have the same file name stub!\n" }

	# read in first record of single-end file
	my $header_s = <$fh_single>;
	my $seq_s = <$fh_single>; 
	print  length($seq_s)-1, "\t";
	my $qual_head_s = <$fh_single>;
	my $qual_s = <$fh_single>;
	if ( !$header_s || !$seq_s || !$qual_head_s || !$qual_s ) { 
		die "$fname_s might be empty or not in right format.\n"; 
	}	
                                                                                                        
	# check the input format of the first record of the single end file
	die "input file not in .fastq format: $fname_s\n" if
	( 
		$header_s !~ /^@.+$/
	  	or $seq_s =~ /^[^AGCTN]$/
	  	or $qual_head_s !~ /^\+/
	  	or $qual_s =~ /^[^!-h]$/
	);
                                                                                                        
	# read in first record of paired-end file
	my $header_p = <$fh_paired>;	
	my $seq_p = <$fh_paired>; 
	print  length($seq_p)-1, "\n";
	my $qual_head_p = <$fh_paired>;
	my $qual_p = <$fh_paired>;
	
	if ( !$header_p || !$seq_p || !$qual_head_p || !$qual_p ) { 
		die "$fname_p might be empty or not in right format.\n"; 
	}

	if ( substr($header_s, 0, -2) ne substr($header_p, 0, -2) ) {
		die "paired files do not have corresponding records in the same order!\n$header_s\t$header_p\n";                     
	}

	# check the input format of the first record of the paired end file
	die "input file not in .fastq format: $fname_p\n" if
       	( 
		$header_p !~ /^@.+$/
         	or $seq_p =~ /^[^AGCTN]$/
         	or $qual_head_p !~ /^\+/
         	or $qual_p =~ /^[^!-h]$/
       	);        
                                                                                                        
	# open output files
	open (my $single_out, ">", "$fname_s_stub\_trueTags.fq_1") 
		or die "Can't write to file $fname_s_stub\_trueTags.fq_1: $!";
	open (my $paired_out, ">", "$fname_p_stub\_trueTags.fq_2")
		or die "Can't write to file $fname_p_stub\_trueTags.fq_2: $!";
	
	# if specified on the command line, 
	# search for only exact matches
	if ($exact) { 
		($discarded, $kept) = filter_exact_matches(
			{
				single_out => $single_out,
				paired_out => $paired_out,
				fh_single => $fh_single,
				fh_paired => $fh_paired,
				header_s => $header_s,
				seq_s => $seq_s,
				qual_head_s => $qual_head_s,
				qual_s => $qual_s,
				header_p => $header_p,
				seq_p => $seq_p,
				qual_head_p => $qual_head_p,
				qual_p => $qual_p,
			}
		);
	# search for exact and close matches
	# of the search string with the beginning
	# of the sequences
	}else { 
		($discarded, $kept) = filter_close_matches(
			{
				single_out => $single_out,
				paired_out => $paired_out,
				fh_single => $fh_single,
				fh_paired => $fh_paired,
				header_s => $header_s,
				seq_s => $seq_s,
				qual_head_s => $qual_head_s,
				qual_s => $qual_s,
				header_p => $header_p,
				seq_p => $seq_p,
				qual_head_p => $qual_head_p,
				qual_p => $qual_p,
			}
		);
	}
	if ($kept) { return $discarded, $kept, $fname_s_stub; }
	else { return $discarded, $discarded, $fname_s_stub; }
}

########################################################################
sub parse_command_line {
	#die $usage unless @ARGV;
	while(@ARGV){
		$_ = shift @ARGV;
		if ( $_ =~ /^-s$/ ) { $glob_single= shift @ARGV; }
		if ( $_ =~ /^-p$/ ) { $glob_paired = shift @ARGV; }
		if ( $_ =~ /^-t$/ ) { $max_processes = shift @ARGV; }
		if ( $_ =~ /^-f1$/ ) { $pattern_s = shift @ARGV; }
		if ( $_ =~ /^-f2$/ ) { $pattern_p = shift @ARGV; }
		if ( $_ =~ /^-e$/ ) { $exact = "True"; }
		if ( $_ =~ /^-m$/ ) { $mismatch = shift @ARGV; }
		if ( $_ =~ /^-pos$/ ) { $start_position = (shift @ARGV)-1; }
		if ( $_ =~ /^-inv$/ ) { $invert_match = "True"; }
		if ( $_ =~ /^-h$/ ) { die $usage; }
	}
}

##########################################################################
# FUNCTION: filter from a set of paired sequence files the records where
#	    the search string is found at the beginning of the sequence.
#	    Allow for the user specified number of mismatches.
# RECEIVES: the two file handles of the paired files, the two file handles
#	    for the output files, the first two records of the paired files 
#	    and references to the counting variables $discarded and $kept
# RETURNS:  nothing, prints to file
##########################################################################
sub filter_close_matches {

	my $arg_ref = shift;

	my $single_out = $arg_ref->{single_out};
	my $paired_out = $arg_ref->{paired_out};
	my $fh_single = $arg_ref->{fh_single};
	my $fh_paired = $arg_ref->{fh_paired};
	my $header_s = $arg_ref->{header_s};
	my $seq_s = $arg_ref->{seq_s};
	my $qual_head_s = $arg_ref->{qual_head_s};
	my $qual_s = $arg_ref->{qual_s};
	my $header_p = $arg_ref->{header_p};
	my $seq_p = $arg_ref->{seq_p};
	my $qual_head_p = $arg_ref->{qual_head_p};
	my $qual_p = $arg_ref->{qual_p};
	my $kmer_s = "";
	my $kmer_p = "";
	my $diff_s = 0;
	my $diff_p = 0;
	my $discarded = 0;
	my $kept = 0;

	if ( $pattern_s && $pattern_p ) {
		$kmer_s = substr($seq_s, $start_position, length($pattern_s));
		$kmer_p = substr($seq_p, $start_position, length($pattern_p));
		$diff_s = ($kmer_s ^ $pattern_s) =~ tr[\001-255][];
		$diff_p = ($kmer_p ^ $pattern_p) =~ tr[\001-255][];
		if ($invert_match) {
			unless ( $diff_s <= $mismatch and $diff_p <= $mismatch ) { 
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}else {
			if ( $diff_s <= $mismatch and $diff_p <= $mismatch ) { 
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}
	}
	elsif( $pattern_s ) {
		$kmer_s = substr($seq_s,  $start_position, length($pattern_s));
		$diff_s = ($kmer_s ^ $pattern_s) =~ tr[\001-\255][];
		if ($invert_match) {
			unless ( $diff_s <= $mismatch ) { 
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}else {
			if ( $diff_s <= $mismatch ) { 
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}
	}
	elsif ( $pattern_p ) {
		$kmer_p = substr($seq_p, $start_position, length($pattern_p));
		$diff_p = ($kmer_p ^ $pattern_p) =~ tr[\001-\255][];
		if ($invert_match) {
			unless ( $diff_p <= $mismatch ) { 
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}else {
			if ( $diff_p <= $mismatch ) { 
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}
	}else{ die "please provide at least one search pattern\n", $usage; }
	

	# now read in the rest of the two paired files
	while ( defined($header_s = <$fh_single>)
			and defined($header_p = <$fh_paired>) )
	{
		$seq_s = <$fh_single>;
		$qual_head_s = <$fh_single>;
		$qual_s = <$fh_single>;
	
		$seq_p = <$fh_paired>;
		$qual_head_p = <$fh_paired>;
		$qual_p = <$fh_paired>;
	
		if ( substr($header_s, 0, -2) ne substr($header_p, 0, -2) ) {
			die "paired files do not have corresponding records in the same order!\n$header_s\t$header_p\n"     	    
		}

		if ( $pattern_s && $pattern_p ) {
			$kmer_s = substr($seq_s, $start_position, length($pattern_s));
			$kmer_p = substr($seq_p, $start_position, length($pattern_p));
			$diff_s = ($kmer_s ^ $pattern_s) =~ tr[\001-255][];
			$diff_p = ($kmer_p ^ $pattern_p) =~ tr[\001-255][];
			if ($invert_match) {
				unless ( $diff_s <= $mismatch and $diff_p <= $mismatch ) { 
					print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
					print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
					$kept++;
				}else { $discarded++; }
			}else {
				if ( $diff_s <= $mismatch and $diff_p <= $mismatch ) { 
					print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
					print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
					$kept++;
				}else { $discarded++; }
			}
#			if ( $diff_s <= $mismatch  and $diff_p <= $mismatch ) { 
#				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
#				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
#				$kept++;
#			}else { $discarded++; }
		}
		elsif( $pattern_s ) {
			$kmer_s = substr($seq_s, $start_position, length($pattern_s));
			$diff_s = ($kmer_s ^ $pattern_s) =~ tr[\001-\255][];
			if ($invert_match) {
				unless ( $diff_s <= $mismatch ) { 
					print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
					print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
					$kept++;
				}else { $discarded++; }
			}else {
				if ( $diff_s <= $mismatch ) { 
					print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
					print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
					$kept++;
				}else { $discarded++; }
			}
			#if ( $diff_s <= $mismatch ) { 
			#	print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
			#	print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
			#	$kept++;
			#}else { $discarded++; }
		}
		elsif ( $pattern_p ) {
			$kmer_p = substr($seq_p, $start_position, length($pattern_p));
			$diff_p = ($kmer_p ^ $pattern_p) =~ tr[\001-\255][];
			if ($invert_match) {
				unless ( $diff_p <= $mismatch ) { 
					print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
					print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
					$kept++;
				}else { $discarded++; }
			}else {
				if ( $diff_p <= $mismatch ) { 
					print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
					print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
					$kept++;
				}else { $discarded++; }
			}
			#if ( $diff_p <= $mismatch ) { 
			#	print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
			#	print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
			#	$kept++;
			#}else { $discarded++; }
		}
	}
	return $discarded, $kept;
}

##############################################################################
# FUNCTION: filter from a set of paired sequence files the records where
#	    the exact search string is found at the beginning of the sequence.
# RECEIVES: the two file handles of the paired files, the two file handles
#	    for the output files, the first two records of the paired files 
#	    and references to the counting variables $discarded and $kept
# RETURNS:  nothing, #prints to file
##############################################################################
sub filter_exact_matches {

	my $arg_ref = shift;

	my $single_out = $arg_ref->{single_out};
	my $paired_out = $arg_ref->{paired_out};
	my $fh_single = $arg_ref->{fh_single};
	my $fh_paired = $arg_ref->{fh_paired};
	my $header_s = $arg_ref->{header_s};
	my $seq_s = $arg_ref->{seq_s};
	my $qual_head_s = $arg_ref->{qual_head_s};
	my $qual_s = $arg_ref->{qual_s};
	my $header_p = $arg_ref->{header_p};
	my $seq_p = $arg_ref->{seq_p};
	my $qual_head_p = $arg_ref->{qual_head_p};
	my $qual_p = $arg_ref->{qual_p};
	my $discarded = 0;
	my $kept = 0;

	# if a search pattern for the single end reads and the paired-end
	# reads has been been provided on the command line
	if ( $pattern_s && $pattern_p ) {
		if ($invert_match) {
			unless ( $seq_s =~ /$pattern_s/o and $seq_p =~ /$pattern_p/o ) { 
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}else {
			if ( $seq_s =~ /$pattern_s/o and $seq_p =~ /$pattern_p/o ) { 
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}
	}
	elsif ( $pattern_s ) {
		if ($invert_match) {
			unless ( $seq_s =~ /$pattern_s/o ) {
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}else {
			if ( $seq_s =~ /$pattern_s/o ) {
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}
	}
	elsif ( $pattern_p ) {
		if ($invert_match) {
			unless ( $seq_p =~ /$pattern_p/o ) {
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}else {
			if ( $seq_p =~ /$pattern_p/o ) {
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}
	}else{ die "please provide at least one search pattern\n", $usage; }
	                                                                                                                                      
	# now read in the rest of the two paired files
	while ( defined($header_s = <$fh_single>)
			and defined($header_p = <$fh_paired>) )
	{
		$seq_s = <$fh_single>;
		$qual_head_s = <$fh_single>;
		$qual_s = <$fh_single>;
	
		$seq_p = <$fh_paired>;
		$qual_head_p = <$fh_paired>;
		$qual_p = <$fh_paired>;
	
		if ( substr($header_s, 0, -2) ne substr($header_p, 0, -2) ) {
			die "paired files do not have corresponding records in the same order!\n$header_s\t$header_p\n"     	    
		}
	                                                                                                                                 
		# if a search pattern for the single end reads and the paired-end
		# reads has been been provided on the command line
		if ( $pattern_s && $pattern_p ) {
			if ($invert_match) {
				unless ( $seq_s =~ /$pattern_s/o and $seq_p =~ /$pattern_p/o ) { 
					print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
					print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
					$kept++;
				}else { $discarded++; }
			}else {
				if ( $seq_s =~ /$pattern_s/o and $seq_p =~ /$pattern_p/o ) { 
					print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
					print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
					$kept++;
				}else { $discarded++; }
			}
		}
		# if only a search pattern for the single-end reads has been provided
		# on the command line
		elsif ( $pattern_s ) {
			if ($invert_match) {
				unless ( $seq_s =~ /$pattern_s/o ) {
					print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
	        			print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
					$kept++;
				}else { $discarded++; }
			}else {
				if ( $seq_s =~ /$pattern_s/o ) {
					print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
	        			print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
					$kept++;
				}else { $discarded++; }
			}
		}
		# if only a search pattern for the paired-end reads has been provided on the 
		# command line
		elsif ( $pattern_p ) {
			if ($invert_match) {
				unless ( $seq_p =~ /$pattern_p/o ) {
					print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
					print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
					$kept++;
				}else { $discarded++; }
			}else {
				if ( $seq_p =~ /$pattern_p/o ) {
					print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
					print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
					$kept++;
				}else { $discarded++; }
			}
		}	
	}
	return $discarded, $kept;
}
