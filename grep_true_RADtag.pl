#!/usr/bin/env perl
use strict; use warnings;

# read in paired file names
# for each pair of .fastq files: read in records, search sequences for remainder of restriction 
# site sequence at the beginning, print out record if true RAD tag sequence

my $usage = "
\nThis script takes a number of paired .fastq file names and removes records if one of their reads
foes not begin with the expected remainder of the restriction site sequence.
The input files must be in fastq format. Quality scores in Phred+64 and Phred+33 are accepted.
Paired files must have corresponding records in the same order and file names must be the same,
except for the file name extension (which should be preceded by the only dot in the file name).

$0

-s single end file name or wildcard expression like \"*.fq_1\"
(wildcards must be protected with quotation marks against shell expansion)
-p paired end file name or wildcard expression like \"*.fq_2\"
-f1 pattern to search for in the single end sequences (can include perl regular expressions)
-f2 pattern to search for in the paired end sequences
-t number of processors available
\n";

# Non-core modules
eval {
     require Parallel::ForkManager;
     Parallel::ForkManager->import();
     };
die "remove_uncalled_bases.pl requires the CPAN module Parallel::ForkManager for parallelising. Please install this package and add it to your Perl library path.\n" if $@;

# initialize variables
my $max_processes = 1;
my ($glob_single, $glob_paired, $pattern_s, $pattern_p, $filtered);
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

# initialize parallel processing
my $pm = new Parallel::ForkManager( scalar(@single_end_files) );
$pm->set_max_procs($max_processes);


# iterate over the paired files
foreach my $file ( 0..$#single_end_files ) {
       print "Scanning reads from files: ", $single_end_files[$file], "\t", $paired_end_files[$file], "\n";
       #print $single_end_files[$file], "\t", $paired_end_files[$file], "\n";

       # start this loop for as many filenames as there are processors available
       $pm->start and next;

       open (my $single_in, "<", $single_end_files[$file]) or die "Can't open " . $single_end_files[$file] . "$!\n";
       open (my $paired_in, "<", $paired_end_files[$file]) or die "Can't open " . $paired_end_files[$file] . "$!\n";

       my ($filtered, $fname_s_stub) = filter_true_tags (
		{
			fh_single => $single_in,
			fh_paired => $paired_in,
			single_file_name => $single_end_files[$file],
			paired_end_name => $paired_end_files[$file]
		}
       );	 
       printf STDERR "%s%.2f%s", "discarded ", $filtered, "% of records from sample \"$fname_s_stub\"\n"; 
       $pm->finish;       
}
$pm->wait_all_children;

#############################################################################
### FUNCTION: searches sequences in paired files for specified patterns
###	      and prints only those records where the pattern has been found 
### RECEIVES: file handles and file names of paired files to process
### RETURNS: percentage discarded records, sample name
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
	unless ( ($fname_s_stub ) = ($fname_s =~ /(.+)\..+/) ) { die $!; }	
	unless ( ($fname_p_stub ) = ($fname_p =~ /(.+)\..+/) ) { die $!; }
	
	if ( $fname_s_stub ne $fname_p_stub ) { die "paired files don't have the same file name stub!\n" }

	# read in first record of single-end file
	my $header_s = <$fh_single>;
	my $seq_s = <$fh_single>; 
	print length($seq_s)-1, "\t";
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
	print length($seq_p)-1, "\n";
	my $qual_head_p = <$fh_paired>;
	my $qual_p = <$fh_paired>;
	
	if ( !$header_p || !$seq_p || !$qual_head_p || !$qual_p ) { 
		die "$fname_p might be empty or not in right format.\n"; 
	}

	if ( substr($header_s, 0, -2) ne substr($header_p, 0, -2) ) {
		die "paired files do not have corresponding records in the same order!\n$header_s\t$header_p\n";                     }

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
	open(my $paired_out, ">", "$fname_p_stub\_trueTags.fq_2")
		or die "Can't write to file $fname_p_stub\_trueTags.fq_2: $!";

	# if a search pattern for the single end reads and the paired-end
	# reads has been been provided on the command line
	if ( $pattern_s && $pattern_p ) {
		if ( $seq_s =~ /$pattern_s/o and $seq_p =~ /$pattern_p/o ) { 
			print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
			print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
			$kept++;
		}else { $discarded++; }
	}
	elsif ( $pattern_s ) {
		if ( $seq_s =~ /$pattern_s/o ) {
			print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
                	print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
			$kept++;
		}else { $discarded++; }
		
	}
	elsif ( $pattern_p ) {
		if ( $seq_p =~ /$pattern_p/o ) {
			print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
			print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
			$kept++;
		}else { $discarded++; }
	}

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
			die "paired files do not have corresponding records in the same order!\n$header_s\t$header_p\n"     		     }
        
		# if a search pattern for the single end reads and the paired-end
		# reads has been been provided on the command line
		if ( $pattern_s && $pattern_p ) {
			if ( $seq_s =~ /$pattern_s/o and $seq_p =~ /$pattern_p/o ) { 
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}
		# if only a search pattern for the single-end reads has been provided
		# on the command line
		elsif ( $pattern_s ) {
			if ( $seq_s =~ /$pattern_s/o ) {
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
                		print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		
		}
		# if only a search pattern for the paired-end reads has been provided on the 
		# command line
		elsif ( $pattern_p ) {
			if ( $seq_p =~ /$pattern_p/o ) {
				print $single_out $header_s, $seq_s, $qual_head_s, $qual_s; 
				print $paired_out $header_p, $seq_p, $qual_head_p, $qual_p;
				$kept++;
			}else { $discarded++; }
		}
	
	}
	return($discarded*100/$kept, $fname_s_stub);
}

########################################################################
sub parse_command_line {
	#die $usage unless @ARGV;
	while(@ARGV){
		$_ = shift @ARGV;
		if ( $_ =~ /^-s/ ) { $glob_single= shift @ARGV; }
		if ( $_ =~ /^-p$/ ) { $glob_paired = shift @ARGV; }
		if ( $_ =~ /^-t$/ ) { $max_processes = shift @ARGV; }
		if ( $_ =~ /^-f1/ ) { $pattern_s = shift @ARGV; }
		if ( $_ =~ /^-f2/ ) { $pattern_p = shift @ARGV; }
	}
}


















































































































































































































































































































































































































