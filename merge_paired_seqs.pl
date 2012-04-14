#!/usr/bin/env perl
# merge_paired_seqs.pl
use strict; use warnings;

# get single and paired end file names
# for each sample open paired sequence files
# subroutine 'merge': read in paired files record by record, 
# print out merged records to new file
#
my $usage = "
This programme takes matching single-end and paired-end read files and merges their sequences.
It creates a new fastq file containing records with concatenated sequences and quality scores.
The input files must be in fastq format. Quality scores in Phred+64 and Phred+33 are accepted.
Paired files must have corresponding records in the same order and file names must be the same,
except for the file name extension (which should be preceded by the only dot in the file name).

$0

-s single-end read file name or wildcard expression like \"*.fq_1\"
(wildcards must be protected with quotation marks against shell expansion)
-p paired-end read file name or wildcard expression like \"*.fq_2\" 
-t the number of threads to run (default: 1)
\n";

# Non-core modules
eval {
     require Parallel::ForkManager;
             Parallel::ForkManager->import();
     };
die "merge_paired_seqs.pl requires the CPAN module Parallel::ForkManager for parallelising. Please install this package and add it to your Perl library path.\n" if $@;

my $max_processes = 1;
my ($glob_single, $glob_paired);
my @arg = @ARGV;

parse_command_line();
unless ( defined($glob_single) && defined($glob_paired) ) { 
	print STDERR "\nyour command line:\n", $0, " ",  "@arg\n";
	print STDERR $usage; 
	exit; 
}

my @single_end_files = glob($glob_single);
my @paired_end_files = glob($glob_paired);

# if no single-end or paired-end files are found 
unless ( @single_end_files and @paired_end_files ) { 
	print STDERR  "\nyour command line:\n", $0, " ",  "@arg\n";
	print STDERR  $usage; 
	exit; 
}

print "command line used for this run:\n", $0, " ", "@arg\n", $usage;

# initialise the parallel processing
my $pm = new Parallel::ForkManager( scalar(@single_end_files) );
$pm->set_max_procs($max_processes);

#  sort the lists of input files
@single_end_files = sort @single_end_files;
@paired_end_files = sort @paired_end_files;

# iterate over the paired files
foreach my $file ( 0..$#single_end_files ) {
	print "Merging reads from files: ", $single_end_files[$file], "\t", $paired_end_files[$file], "\n";
	#print $single_end_files[$file], "\t", $paired_end_files[$file], "\n";

	# start this loop for as many filenames as there are processors available
	$pm->start and next;

	open (my $single_in, "<", $single_end_files[$file]) or die "Can't open " . $single_end_files[$file] . "$!\n";
	open (my $paired_in, "<", $paired_end_files[$file]) or die "Can't open " . $paired_end_files[$file] . "$!\n";

	merge (
		{
		single_fh => $single_in,
		paired_fh => $paired_in,
		single_file_name => $single_end_files[$file],
		paired_end_name => $paired_end_files[$file]
		}

	);
	$pm->finish;
}
$pm->wait_all_children;
print STDERR "Finished\n";

#################################################################
### FUNCTION: merge single-end and paired-end read of paired files
###	      and print them out to a new file
### RECEIVES: filehandles and filenames of paired files 
### RETURNS:  nothing, prints to file
#################################################################
sub merge {
	my $arg_ref = shift;

	# initialize variables
	my $fh_single = $arg_ref->{single_fh};
	my $fh_paired = $arg_ref->{paired_fh};
	my $fname_s = $arg_ref->{single_file_name};
	my $fname_p = $arg_ref->{paired_end_name};
	my ($fname_s_stub, $fname_p_stub);


	# get sample names
 	unless ( ($fname_s_stub ) = ($fname_s =~ /(.+)\..+/) ) {
		die $!;
	}
	unless ( ($fname_p_stub ) = ($fname_p =~ /(.+)\..+/) ) {
		die $!;
	}
	if ( $fname_s_stub ne $fname_p_stub ) { die "paired files don't have the same file name stub!\n" }

	# read in first record of single-end file
	my $header_s = <$fh_single>;
	$header_s = substr $header_s, 0, -2;
	chomp( my $seq_s = <$fh_single> ); print length($seq_s), "\t";
	my $qual_head_s = <$fh_single>;
	chomp( my $qual_s = <$fh_single> );
	if ( !$header_s || !$seq_s || !$qual_head_s || !$qual_s ) { 
		die "$fname_s might be empty or not in right format.\n"; 
	}

	# check the input format of the first record of the single end file
	die "input file not in .fastq format: $fname_s 
		$header_s\n" if
		( $header_s !~ /^@.+$/
		  or $seq_s =~ /^[^AGCTN]$/
		  or $qual_head_s !~ /^\+/
		  or $qual_s =~ /^[^!-h]$/
		);

	# read in first record of paired-end file
	my $header_p = <$fh_paired>;
	$header_p = substr $header_p, 0, -2;
	my $seq_p = <$fh_paired>; print length($seq_p)-1, "\n";
	my $qual_head_p = <$fh_paired>;
	my $qual_p = <$fh_paired>;
	if ( !$header_p || !$seq_p || !$qual_head_p || !$qual_p ) { 
		die "$fname_p might be empty or not in right format.\n"; 
	}
	if ( $header_s ne $header_p ) {
		die "paired files do not have corresponding records in the same order!\n$header_s\t$header_p\n"

	}

	# check the input format of the first record of the paired end file
	die "input file not in .fastq format: $fname_p\n" if
	       ( $header_p !~ /^@.+$/
	         or $seq_p =~ /^[^AGCTN]$/
	         or $qual_head_p !~ /^\+/
	         or $qual_p =~ /^[^!-h]$/
	       );        

	# open output files
	open (my $merged_out, ">", "$fname_s_stub\_merged.fq") 
		or die "Can't write to file $fname_s_stub\_merged.fq: $!";

	# merge and print out the first paired record
	print $merged_out
		$header_s . "merged\n",
		$seq_s . $seq_p,
		$qual_head_s,
		$qual_s . $qual_p;

	# read in the rest of the two files and merge paired records
	while ( defined($header_s = <$fh_single>) 
		and defined($header_p = <$fh_paired>) )
	{

		$header_s = substr $header_s, 0, -2;
		chomp( $seq_s = <$fh_single> );
		$qual_head_s = <$fh_single>;
		chomp( $qual_s = <$fh_single> );
        	die "$header_s\n$seq_s$qual_head_s$qual_s\n$!" 
			if ( !$header_s || !$seq_s || !$qual_head_s || !$qual_s );
                                                                              
		$header_p = substr $header_p, 0, -2;
		$seq_p = <$fh_paired>;
		$qual_head_p = <$fh_paired>;
		$qual_p = <$fh_paired>;
		if ( !$header_p || !$seq_p || !$qual_head_p || !$qual_p ) {
			die "paired-end file contains fewer records than the single-end file
			of sample: $fname_s_stub", "\n";
		}
		# check whether paired files have corresponding records
		if ( $header_s ne $header_p ) {
			die "paired files do not have corresponding records in the same order!\nLine: $.\n$header_s\t$header_p\n";
		}			

		# merge paired records and print out to file
		print $merged_out
			$header_s . "merged\n",
			$seq_s . $seq_p,
			$qual_head_s,
			$qual_s . $qual_p;
	}
}

#################################################################

sub parse_command_line {
	die $usage unless @ARGV;
	while(@ARGV){
		$_ = shift @ARGV;
		if ( $_ =~ /^-p$/ ) { $glob_paired = shift @ARGV; }
		if ( $_ =~ /^-s$/ ) { $glob_single = shift @ARGV; }
		if ( $_ =~ /^-t$/ ) { $max_processes = shift @ARGV; }
		if ( $_ =~ /^-h$/ ) { die $usage; }
	}
}

































































































































































































