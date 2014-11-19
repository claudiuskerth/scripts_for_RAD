#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: find_linked_RADtags.pl
#
#        USAGE: ./find_linked_RADtags.pl  
#
#  DESCRIPTION: This script is designed to search a sorted bam file created from 
#  				standard RAD paired-end read data for reference contigs that have read pairs 
#  				mapped to both sides of a single SbfI restriction site (or any other
#  				restriction enzyme cut site creating a 4 bp overhang).
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 23/01/14 16:19:54
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
#use diagnostics;

my @which = `which samtools`;
die "samtools needs to be installed and in your PATH\n" unless @which;

my $usage = "
$0 [options] <input1.bam> <input2.bam> ...

This script is designed to search a sorted bam file created by mapping 
standard RAD paired-end read data for reference contigs that have read pairs 
mapped to both sides of a single SbfI restriction site (or any other
restriction enzyme cut site creating an overhang).
Be aware that this script expects input files that were position sorted with the
\'samtools sort\' command.

-[o|overlap] [INT]  --> sticky end length that restriction enzyme leaves (default: 4)
-[p|paired_end]     --> collect PE reads upstream and downstream of detected RAD sites into separate files 
-[h|help]           --> print this help message
-[v|verbose]        --> print more verbose run reports

In order to detect a potential linked RADtags site, this script requires at least one input file 
to provide one properly mapped read pair on each side of a potential restriction site with their 
SE reads overlapping by the specified length. Only the most 5 prime RADtag site in a contig 
is currently detected.
\n";

die $usage unless @ARGV;

my $overlap = 4;
my @inputFiles;
my $verbose = 0;
my $collect_PE_reads = 0;

parse_command_line();

if($verbose){
	print "Setting restriction site overhang to $overlap.\n";
	print "The following input files will be scanned:\n";
	foreach my $input (@inputFiles){ print "$input\n"; }
}

# global variables
my ($flag, %map_pos, @fields, $contig_ID, $cigar, $map_pos, $stub, $header, $pp);

for my $filename (@inputFiles){ 

	#-------------------------------------------------------------------------------
	#  file tests
	#-------------------------------------------------------------------------------
	die "File $filename does not exist\n" unless (-e $filename);
	die "$filename is a directory\n" if (-d $filename);
	die "File $filename is empty\n" if (-z $filename);

	#-------------------------------------------------------------------------------
	#  open the file
	#-------------------------------------------------------------------------------
	my $IN = pipe_open($filename);

	#------------------------------------------
	# find contigs with linked RAD tags
	#------------------------------------------
	my $linked_RADtag_contigs_ref = find_contigs_with_linked_RADtags($IN);

	close $IN;

	#-------------------------------------------------------------------------------
	# print the detected contigs to a report file
	#-------------------------------------------------------------------------------
	log_detected_contigs($filename, $linked_RADtag_contigs_ref);

	# index the input file for retrieval of reads with samtools
	#system("samtools index $filename")==0 or die $!;
	#system("samtools index $filename")==0 or die $!;
	# unfortunately, in the collect PE reads step, I also need to collect discordantly mapped PE reads
	# whose SE read mapped to a detected RADtag site. So 'samtools view infile UNIGENE' will
	# be faster but will not get all reads I want.

	if($collect_PE_reads){
		#----------------------------------------------------
		# get paired-end reads from detected RADtag sites
		#----------------------------------------------------
		# foreach contig that contains a detected RADtag site
		foreach my $detected_contig (keys %$linked_RADtag_contigs_ref){
			# get all reads that map to this contig. Note, that unmapped reads
			# get the same mapping position as their mapped mates.
			# The awk command insures that this also gets reads that mapped to other
			# contigs, but whose mate mapped to this contig.
			open( $IN, "samtools view $filename | awk \'/$detected_contig/\' | ") or die $!;

			if($verbose){ print "Reading in bam records for reference contig ", $detected_contig, "\n"; }

			#--------------------------------------------------------------
			#  collect all SE reads that map to the detected RADtag site 
			#  as well as their PE reads
			#--------------------------------------------------------------
			my ($SE_ref, $PE_ref) = collect_all_reads($IN, $detected_contig, $linked_RADtag_contigs_ref);
			
			#---------------------------------------------
			# print out all SE and PE reads that map to a 
			# RAD tag site
			#---------------------------------------------
			print_out_reads($detected_contig, $SE_ref, $PE_ref);

		}# END foreach detected contig
	}
}# END foreach input file

#---------------------------------------------------------------------------

#------------------------
# SUBROUTINES
#------------------------

sub parse_command_line{
	while(@ARGV){
		$_ = shift @ARGV;
		if(/^-+h$/){croak $usage;}
		elsif(/^-+o/){ $overlap = shift @ARGV; }
		elsif(/^-+v$/){ $verbose = 1; }
		elsif(/^-+p/){ $collect_PE_reads = 1; }
		else{ push @inputFiles, $_ ; }
	}
}


#===  FUNCTION  ================================================================
#         NAME: PIPE_OPEN
#      PURPOSE: open a file with an external command an read from its STDOUT
#   PARAMETERS: ????
#  DESCRIPTION: an indirect filehandle for reading
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub pipe_open {
	my $filename = shift;
	my $in;
	if($filename =~ /bam$/){
		# the "-F128" to samtools lets it output only single end reads
		# the awk command takes only read pairs with the "proper pair" flag set
		# i. e. ONLY PROPERLY MAPPED READ PAIRS ARE USED FOR CONTIG DETECTION.
		open( $in, "samtools view -F128 $filename | awk \'and(\$2, 2)\' |" ) 
			or die "Could not fork samtools: $!"; 
#	}elsif($filename =~ /\.sam\./){
#		#open( $IN, "samtools view -S $filename | head -n 10000 | " ) or die $!; 
#		#open( $IN, "samtools view -SF128 $filename | awk \'and(\$2, 2)\' | " ) or die $!; 
#		die "Your filename file needs to end either with bam or sam or sam.gz.\n";
	}else{
		die "Your filename file needs to end with bam and be a position sorted BAM file.\n";
	}
	return $in;
}

#===  FUNCTION  ================================================================
#         NAME: COMP_MAP_POS
#      PURPOSE: Checks for whether any possible pair of forward mapping and reverse
#               mapping reads overlap by the specified length. If yes, then the
#               name of the reference contig is store in a hash together with mapping
#               position of the forward mapping read.
#   PARAMETERS: reference to hash that stores names of contigs that contain a detected
#               linked RADtags site 
#      RETURNS: notihin
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub comp_map_pos{
	my $linked_RADtag_contigs_ref = shift;

	MAPPOS: 
	foreach my $pos1 (@{ $map_pos{F} }){        # %map_pos is global, but re-initialised for each new contig
		foreach my $pos2 (@{ $map_pos{R} }){
			if(($pos2 - $overlap + 1) == $pos1){ 
				$linked_RADtag_contigs_ref->{$contig_ID} = $pos1;
				last MAPPOS; # Store the reference contig name and position and stop searching as soon 
				# as a linked RADtags site is found. Since the input bam is position sorted, only the 
				# most upstream linked RADtags site will thus be considered.
			}
		}
	}	
}


#===  FUNCTION  ================================================================
#         NAME: GET_3PRIME_MAP_POS
#      PURPOSE: gets 3 prime end mapping position of a read
#   PARAMETERS: ????
#      RETURNS: ????
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub get_3prime_map_pos {
	my @shift = $cigar =~ /(\d+)[MI]/g; # taking numbers in CIGAR string makes it independent of read length
	# add sum of length of M and I to mapping position
	map { $map_pos += $_ } @shift;
	# if the read starts with an insertion, the first matching base
	# determines the mapping position of the read. This determines
	# the mapping position of the first base of the read as if the read were
	# mapped without the insertion
	if ($cigar =~ /^(\d+)I/){ $map_pos -= $1; }
	$map_pos--;
}


#===  FUNCTION  ================================================================
#         NAME: FIND_CONTIGS_WITH_LINKED_RADTAGS
#      PURPOSE: 
#   PARAMETERS: indirect filehandle
#      RETURNS: ????
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub find_contigs_with_linked_RADtags {

	my $IN = shift;

	# initialize hash that stores contigs with linked RAD tags
	my %linked_RADtag_contigs;

	# get the first SAM record
	my $new_SAM_record = <$IN>;
	
	# split the SAM fields into an array
	@fields = split('\t', $new_SAM_record);                                            # why shouldn't these variables be private?
	$flag = $fields[1];                                                                # 
	$cigar = $fields[5];                                                               # 
	$map_pos = $fields[3];                                                             # 


	CONTIG: while($new_SAM_record){ # while the SAM record is defined
		# for each new contig initialise a hash storing mapping positions of
		# forward mapping and reverse mapping reads separately
		%map_pos = ();                          # global variable
		# store contig ID
		$contig_ID = $fields[2]; 

		store_mapping_pos();

		# while new SAM records can be read in
		while($new_SAM_record = <$IN>){
			# split the record into fields
			@fields = split('\t', $new_SAM_record);
			$flag = $fields[1];
			$cigar = $fields[5];
			$map_pos = $fields[3];

			# if the SAM record belongs to the next contig
			if($fields[2] ne $contig_ID){
				comp_map_pos(\%linked_RADtag_contigs);
				next CONTIG;
			}else{
				store_mapping_pos();
			}

		}# -- END of File
		# the end of the input file has been reached, thus also the end
		# of records for the last contig. So one last time, check for linked
		comp_map_pos(\%linked_RADtag_contigs);
	}# -- END of CONTIG
	return \%linked_RADtag_contigs;
}


#===  FUNCTION  ================================================================
#         NAME: STORE_MAPPING_POS
#      PURPOSE: store a mapping position for forward and backward reads separately
#   PARAMETERS: changes global variables
#      RETURNS: ????
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub store_mapping_pos {
	# if the reads maps to the reverse strand
	# and does not end with an insertion
	if($flag & 16){
		# the new mapping position is the most 3 prime reference position 
		# that is covered by the read 
		get_3prime_map_pos();
		# store the mapping position under reverse reads
		push @{ $map_pos{R} }, $map_pos; 
	# the read maps to the forward strand
	# and if the read mapping does not start with an insertion
	}else{
		# if the read starts with an insertion, the first matching base
		# determines the mapping position of the read. This determines
		# the mapping position of the first base of the read as if the read were
		# mapped without the insertion:
		if ($cigar =~ /^(\d+)I/){ $map_pos -= $1; }
		push @{ $map_pos{F} }, $map_pos; 
	}
}


#===  FUNCTION  ================================================================
#         NAME: log_detected_contigs
#      PURPOSE: write contig IDs with detected RADtag sites to file
#   PARAMETERS: ????
#      RETURNS: ????
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub log_detected_contigs {
	my	($filename, $linked_RADtag_contigs_ref)	= @_;

	# if the filename is somehwere else than the current directory,
	# then this removes removes file path from the file name given.
	# output files will be in the current working directory
	($stub = $filename) =~ s/.*\/(.*)\.bam/$1/; 

	open( my $OUT, ">", $stub . "_linked_RADtag_contigs.out" ) or die $!;

	print $OUT "Contig", "\t", "Position", "\n";
	foreach my $contig (keys %$linked_RADtag_contigs_ref){
		print $OUT $contig, "\t", $linked_RADtag_contigs_ref->{$contig}, "\n";
	}
	close $OUT;

	if($verbose){ print "Found ", scalar keys %$linked_RADtag_contigs_ref, " reference contigs with linked RADtags\n"; }
	return ;
} ## --- end sub log_detected_contigs



#===  FUNCTION  ================================================================
#         NAME: collect_all_reads
#      PURPOSE: collect all reads from detected RADtag site
#   PARAMETERS: ????
#      RETURNS: ????
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub collect_all_reads {
	my	($IN, $detected_contig, $linked_RADtag_contigs_ref)	= @_;
	# initialise variables
	my ($read_ID, $seq, $qual, %PE_init, %SE, %PE, $pp);

	# while reading in all bam records for this contig
	while(<$IN>){
		@fields = split('\t', $_);
		$read_ID = $fields[0];
		$flag = $fields[1];
		$contig_ID = $fields[2];
		$map_pos = $fields[3];
		$cigar = $fields[5];
		$seq = $fields[9];
		$qual = $fields[10];

		# Note, that in SAM files created by stampy, paired reads have exactly the same 
		# read ID in the first column (any 2 at the end of the ID that was there before
		# is replaced by a 1)

		# if read is second in pair, i. e. a PE read
		if($flag & 128){ 						
			# if the read was reverse complemented by the mapping programme
			if($flag & 16){
				# get the read in its original orientation
				($seq = reverse $seq) =~ tr/ACGT/TGCA/;
				$qual = reverse $qual;
			}
			$PE_init{$read_ID} = [$seq, $qual]; # collect all PE reads
		# if read is first in pair, i. e. a SE read
		}else{ 
			# If only the mate of this single-end read mapped to the current contig,
			# skip recording this read. It cannot come from a RAD site that was detected
			# with the above method.
			next if ($contig_ID !~ $detected_contig);
			# if read is printed as reverse complement of input sequence in the SAM input
			if($flag & 16){ 
				# ignore read if it ends in an insertion
				next if ($cigar =~ /\d+I$/);
				# get most 3 prime mapping position
				get_3prime_map_pos();
				($map_pos -= $overlap) += 1;
				($seq = reverse $seq) =~ tr/ACGT/TGCA/;
				$qual = reverse $qual;
				# if the read maps to the linked RADtags site detected for this contig
				if($map_pos == $linked_RADtag_contigs_ref->{$detected_contig}){
					$pp = $flag & 2 ? 1 : 0;
#					if($flag & 2){ $pp = 1 }else{ $pp = 0 }
					$SE{upstream}->{$read_ID} = [$seq, $qual, $pp];
				}
			# if read is mapped on the same strand as the reference sequence
			}else{
				# ignore read if it begins with an insertion
				next if ($cigar =~ /^\d+I/);
				# if the read maps to the linked RADtags site detected for this contig
				if($map_pos == $linked_RADtag_contigs_ref->{$detected_contig}){
					$pp = $flag & 2 ? 1 : 0;
#					if($flag & 2){ $pp = 1 }else{ $pp = 0 }
					$SE{downstream}->{$read_ID} = [$seq, $qual, $pp];
				}
			}
		}
	}# END of reading in bam records for the contig
	
	#----------------------------------------------------
	# only keep PE reads whose mate maps to a RADtag site
	#----------------------------------------------------
#	foreach my $read_ID (keys %$PE_init_ref){
#		# if there exists a SE read with the same read ID
#		if(exists $SE_ref->{upstream}->{$read_ID} or exists $SE_ref->{downstream}->{$read_ID}){
#			$PE{$read_ID} = $PE_init_ref->{$read_ID};
#		}
#	}
			
	# this code, compared to the above version, should save many unnecessary hash lookups
	# since only PE reads from SE reads that map to a detected RADtag site are to be safed.
	foreach my $read_ID (keys %{ $SE{upstream} }){
		$PE{$read_ID} = $PE_init{$read_ID};
	}
	foreach my $read_ID (keys %{ $SE{downstream} }){
		$PE{$read_ID} = $PE_init{$read_ID};
	}

	return (\%SE, \%PE);
} ## --- end sub collect_all_reads


#===  FUNCTION  ================================================================
#         NAME: print_out_reads
#      PURPOSE: 
#   PARAMETERS: ????
#      RETURNS: ????
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub print_out_reads {
	my ( $detected_contig, $SE_ref, $PE_ref ) = @_;
	
	my $out_name = $stub . "_$detected_contig" . "_downstream" . ".fq";

	# first all reads downstream of RAD site in one file
	open( my $FOR, ">", $out_name) or die $!;

	foreach my $read_ID (keys %{ $SE_ref->{downstream} }){
		$pp = "";
		($header = $read_ID) =~ s[1\/1$][]; 
		# insert "pp_" into fastq header if read pair is "proper"
		if($SE_ref->{downstream}->{$read_ID}->[2]){ $pp = "pp_" }
		# print PE read fastq record
		print $FOR "@", $header, $pp, "2\n",
					$PE_ref->{$read_ID}->[0], "\n",
					"+", "\n",
					$PE_ref->{$read_ID}->[1], "\n";
		# print SE read fastq record
		print $FOR "@", $header, "1\n",
					$SE_ref->{downstream}->{$read_ID}->[0], "\n",
					"+", "\n",
					$SE_ref->{downstream}->{$read_ID}->[1], "\n";
	}
	close $FOR;

	$out_name = $stub . "_$detected_contig" . "_upstream" . ".fq";
	# now all reads upstream of RAD site
	open( my $REV, ">", $out_name) or die $!;

	foreach my $read_ID (keys %{ $SE_ref->{upstream} }){
		$pp = "";
		($header = $read_ID) =~ s[1\/1$][]; 
		if($SE_ref->{upstream}->{$read_ID}->[2]){ $pp = "pp_" }
		print $REV "@", $header, $pp, "2\n",
					$PE_ref->{$read_ID}->[0], "\n",
					"+", "\n",
					$PE_ref->{$read_ID}->[1], "\n";
		print $REV "@", $header, "1\n",
					$SE_ref->{upstream}->{$read_ID}->[0], "\n",
					"+", "\n",
					$SE_ref->{upstream}->{$read_ID}->[1], "\n";
	}
}
