#!/usr/bin/env perl
# trim.pl
use strict; use warnings;

# expects fastq format, and is not specififc for a certain quality score format
# makes offset positive but can be set by the user to any length
# if end trim length is not specified, then everything is kept
# makes end trim length provided by the user negative so that it clips that many bases off the end of the sequence
 
my $usage = "
$0
 
-t the number of threads to run (default: 1)
-i file name or wildcard expression like \"*.fq_2\" (wildcards must be protected with quotation marks against shell expansion) 
-b number of bases to trim from beginning of sequence (default 0) 
-e number of bases to trim from end of sequence (default: 0)

expects fastq format, and is not specififc for a certain quality score format
\n";

# Non-core modules
eval {
    require Parallel::ForkManager;
    Parallel::ForkManager->import();
};
die
"trim.pl requires the CPAN module Parallel::ForkManager for parallelising. Please install this package and add it to your Perl library path.\n"
  if $@;

# initialise variables
my $max_processes = 1;
my $glob = ""; # gets pattern of file names of files to process
my $offset = 0;
my $trim_length;

# get command line parameters
parse_command_line();

# get file names to process
my @files_to_trim = glob($glob);
# ignore already trimmed files
@files_to_trim = grep { !/trimmed/ } @files_to_trim;

# initialise the parallel processing
my $pm = new Parallel::ForkManager( scalar(@files_to_trim) );
$pm->set_max_procs($max_processes);

foreach my $file (@files_to_trim) {
	# start this loop for as many filenames as there are processors available
	$pm->start and next;
#	print $file, "\n";

	my ($sample_name, $ext) = $file =~ /^(.+)\.(.+)$/;
	
	# open single end file 
        open my $single_in, '<', "./$file"
          or die "Can't open input file $file!\n";
    
    	trim(
		{
		fh => $single_in, 
		sample_name => $sample_name, 
		file_ext => $ext, 
		trim_offset => $offset, 
		trim_length => $trim_length
		}
	);
        
    $pm->finish;
	
}
$pm->wait_all_children;

############################################

sub trim {

	my $arg_ref = shift;

	my $infile = $arg_ref->{fh};
	my $sample_name = $arg_ref->{sample_name};
	my $offset = $arg_ref->{trim_offset};
	my $length = $arg_ref->{trim_length};
	my $ext = $arg_ref->{file_ext};
	
	# read in the first single end record and check format
    	my $header    = <$infile>;
	my $seq       = <$infile>;
	chomp $seq;
	# unless the trim length is specified by the user, set it to the length of the sequence
	unless ( $length ) { $length = length($seq); }
	die "The argument to -b must not be longer or equal to the sequence length.\n" if ( $offset >= length($seq) );
	$seq 	      = substr $seq, $offset, $length;
	my $qual_head = <$infile>;
	my $qual      = <$infile>;
	chomp $qual;
	$qual	      = substr $qual, $offset, $length;
	die "input file not in .fastq format: $sample_name\n" if 
		(
		$header !~ /^@.+$/ 
		or $seq =~ /^[^AGCTN]$/
		or $qual_head !~ /^\+/
		or $qual =~ /^[^!-h]$/
		);
	return if (!$header) || (!$seq) || (!$qual_head) || (!$qual);
	
	# open output file
	open( my $trimmed, ">", "$sample_name\_trimmed\.$ext") or die "Can't write to file: $sample_name\_trimmed\.$ext, $!";

	print $trimmed 
		$header, 
		$seq, "\n", 
		$qual_head, 
		$qual, "\n";
	
	# read in the rest of the single end records
    	while (<$infile>) {
		$header    = $_; 
		$seq       = <$infile>;
		chomp $seq;
		$seq 	   = substr $seq, $offset, $length;
		$qual_head = <$infile>;
		$qual      = <$infile>;
		chomp $qual;
		$qual	   = substr $qual, $offset, $length;
		last if (!$header) || (!$seq) || (!$qual_head) || (!$qual);
		print $trimmed 
			$header, 
			$seq, "\n", 
			$qual_head, 
			$qual, "\n";
	}

} 

####################################################

sub parse_command_line {
	die $usage unless @ARGV;
	while (@ARGV) {
		$_ = shift @ARGV;
        	if    ($_ =~ /^-t$/) { $max_processes = shift @ARGV; }
        	elsif ($_ =~ /^-i$/) { $glob = shift @ARGV; }
        	elsif ($_ =~ /^-b$/) { $offset = abs(shift @ARGV); }
        	elsif ($_ =~ /^-e$/) { $trim_length  = 0 - abs(shift @ARGV); } 
        	elsif ($_ =~ /^-h$/) { print $usage; }
	}

}
