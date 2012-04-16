#!/usr/bin/env perl
#remove_uncalled_bases.pl
use strict; use warnings;

# get single end and paired end file names
# for each sample open single end and paired end file
# check for uncalled bases in either single or paired end file

my $usage = "
This programme takes single-end and paired-end read files and removes records if either of their reads
contain an uncalled base, i. e. in the output neither single nor paired-end record should contain an uncalled base.
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
die "remove_uncalled_bases.pl requires the CPAN module Parallel::ForkManager for parallelising. Please install this package and add it to your Perl library path.\n" if $@;

my $max_processes = 1;
my ($glob_single, $glob_paired);

parse_command_line();
unless( defined($glob_single) && defined($glob_paired) ) { print STDERR $usage; exit; }

my @single_end_files = glob($glob_single);
my @paired_end_files = glob($glob_paired);

# if no single-end or paired-end files are found 
unless ( @single_end_files and @paired_end_files ) { 
	print STDERR  "\nyour command line:\n", $0, " ",  "@arg\n";
	print STDERR  $usage; 
	exit; 
}
print $0, " ", "@arg\n";

# initialise the parallel processing
my $pm = new Parallel::ForkManager( scalar(@single_end_files) );
$pm->set_max_procs($max_processes);

# sort the lists of input files
@single_end_files = sort @single_end_files;
@paired_end_files = sort @paired_end_files;

foreach	my $file ( 0..$#single_end_files ) {
	print "Processing files: ", $single_end_files[$file], "\t", $paired_end_files[$file], "\n";
	
	#my ($fname_s_stub, $ext_s) = ($single_end_files[$file] =~ m!/*(.+)\.(.+)$!);
	#my ($fname_p_stub, $ext_p) = ($paired_end_files[$file] =~ m!/*(.+)\.(.+)$!);
	#print $fname_s_stub, "\t", $fname_p_stub, "\n";
	#next;
        
	# start this loop for as many filenames as there are processors available
	$pm->start and next;


	open (my $single_in, "<", $single_end_files[$file]) or die "Can't open " . $single_end_files[$file] . "$!\n";
	open (my $paired_in, "<", $paired_end_files[$file]) or die "Can't open " . $paired_end_files[$file] . "$!\n";

	my ($discarded, $kept) = rm_records_with_N (
		{
			single_fh => $single_in,
			paired_fh => $paired_in,
			single_file_name => $single_end_files[$file],
			paired_file_name => $paired_end_files[$file],
		}

	);
	if ($discarded =~ /\D/) { print STDERR $discarded; }
	else {
		printf "%s%d%s%s%s%s%.2f%s", 
			"Discarded ", $discarded, " records from files ", 
			$single_end_files[$file], " and ", "$paired_end_files[$file] (", $discarded*100/$kept, "%)\n";
	}
	$pm->finish;
}
$pm->wait_all_children;
print STDERR "Finished\n";

##################################################################################

sub rm_records_with_N {
	my $arg_ref = shift;
	
	# initialize variables
	my $fh_single = $arg_ref->{single_fh};
	my $fh_paired = $arg_ref->{paired_fh};
	my $fname_s = $arg_ref->{single_file_name};
	my $fname_p = $arg_ref->{paired_file_name};
	my $N_count = 0;
	my ($fname_s_stub, $fname_p_stub);
	my ($ext_s, $ext_p);
	my $discarded = 0;
	my $kept = 0;

	# get sample names
 	unless ( ($fname_s_stub, $ext_s) = ($fname_s =~ m!/*(.+)\.(.+)!) ) {
		die $!;
	}
	unless ( ($fname_p_stub, $ext_p) = ($fname_p =~ m!/*(.+)\.(.+)!) ) {
		die $!;
	}
	if ( $fname_s_stub ne $fname_p_stub ) { die "paired files don't have the same file name stub!\n" }
	
	# read in first record of single-end file
	my $header_s = <$fh_single>;
	my $seq_s = <$fh_single>;
	my $qual_head_s = <$fh_single>;
	my $qual_s = <$fh_single>;
	if ( !$header_s || !$seq_s || !$qual_head_s || !$qual_s ) { return "$fname_s might be empty.\n"; }

	# check the input format of the first record of the single end file
	die "input file not in .fastq format: $fname_s\n" if
		( $header_s !~ /^@.+$/
		  or $seq_s =~ /^[^AGCTN]$/
		  or $qual_head_s !~ /^\+/
		  or $qual_s =~ /^[^!-h]$/
		);

	# read in first record of paired-end file
	my $header_p = <$fh_paired>;
	my $seq_p = <$fh_paired>;
	my $qual_head_p = <$fh_paired>;
	my $qual_p = <$fh_paired>;
	if ( !$header_p || !$seq_p || !$qual_head_p || !$qual_p ) { return "$fname_p might be empty.\n"; }

	if ( substr($header_s, 0, -2) ne substr($header_p, 0, -2) ) {
		die "paired files do not have corresponding records in the same order!\n$header_s\t$header_p\n";
	}

	# check the input format of the first record of the paired end file
	die "input file not in .fastq format: $fname_p\n" if
	       ( $header_p !~ /^@.+$/
	         or $seq_p =~ /^[^AGCTN]$/
	         or $qual_head_p !~ /^\+/
	         or $qual_p =~ /^[^!-h]$/
	       );        

	# open output files
	open (my $single_out, ">", "$fname_s_stub\_cleaned.$ext_s") 
		or die "Can't write to file $fname_s_stub\_cleaned.$ext_s: $!";
	open (my $paired_out, ">", "$fname_p_stub\_cleaned.$ext_p") 
		or die "Can't write to file $fname_p_stub\_cleaned.$ext_p: $!";

	# check whether there are uncalled bases in either 
	# single end or paired end sequence
	$N_count = ($seq_s =~ tr/N/N/) + ($seq_p =~ tr/N/N/);
	
	if ($N_count == 0) {
		print $single_out 
			$header_s, $seq_s, $qual_head_s, $qual_s;
		print $paired_out
			$header_p, $seq_p, $qual_head_p, $qual_p;
		$kept++;
                                                                  
	}else{ $discarded++; }	

	# read in the rest of records in the two paired files in parallel
	while (<$fh_single>){

		$N_count = 0;

		$header_s = $_;
		$seq_s = <$fh_single>;
		$qual_head_s = <$fh_single>;
		$qual_s = <$fh_single>;
	        return if (!$header_s) || (!$seq_s) || (!$qual_head_s) || (!$qual_s);

		$header_p = <$fh_paired>;
		$seq_p = <$fh_paired>;
		$qual_head_p = <$fh_paired>;
		$qual_p = <$fh_paired>;
		return if (!$header_p) || (!$seq_p) || (!$qual_head_p) || (!$qual_p);
		
		$N_count = ($seq_s =~ tr/N/N/) + ($seq_p =~ tr/N/N/);
		
		if ($N_count == 0) {
			print $single_out 
				$header_s, $seq_s, $qual_head_s, $qual_s;
			print $paired_out
				$header_p, $seq_p, $qual_head_p, $qual_p;
			$kept++;

		}else{ $discarded++; }	

	}
	return($discarded, $kept);
}

####################################################################################

sub parse_command_line {
	die $usage unless @ARGV;
	while(@ARGV){
		$_ = shift @ARGV;
		if ( $_ =~ /^-p$/ ) { $glob_paired = shift @ARGV; }
		if ( $_ =~ /^-s$/ ) { $glob_single = shift @ARGV; }
		if ( $_ =~ /^-t$/ ) { $max_processes = shift @ARGV; }
		elsif ($_ =~ /^-h$/ ) { print STDERR $usage; exit; }	
	}
}
