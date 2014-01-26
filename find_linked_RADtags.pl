#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: find_linked_RADtags.pl
#
#        USAGE: ./find_linked_RADtags.pl  
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
#      CREATED: 23/01/14 16:19:54
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

# this programme assumes a sorted sam or bam file

# TODO: add switch to input length of single-end reads

my @which = `which samtools`;
die "samtools needs to be installed and in your PATH\n" unless @which;


my $IN;

for my $input (@ARGV){ # for each input file given
	if($input =~ /bam$/){
		# the "-F128" to samtools lets it output only single end reads
		# the awk command takes only read pairs with the "proper pair" flag set
		open( $IN, "samtools view -F128 $input | awk \'and(\$2, 2)\' | " ) or die $!;
	}elsif($input =~ /\.sam\./){
		#open( $IN, "samtools view -S $input | head -n 10000 | " ) or die $!; 
		open( $IN, "samtools view -SF128 $input | awk \'and(\$2, 2)\' | " ) or die $!; 
	}else{
		die "Your input file needs to end either with bam or sam or sam.gz.\n"
	}

	my ($flag, %map_pos, @fields, $contig_ID, $cigar, $map_pos, @linked_RADtag_contigs);

	# get the first SAM record
	my $new_SAM_record = <$IN>;
	# split the SAM fields into an array
	@fields = split('\t', $new_SAM_record);

	CONTIG:
	while($new_SAM_record){ # while the SAM record is defined
		# initialise a hash storing mapping positions of
		# forward mapping and reverse mapping reads separately
		%map_pos = ();
		$contig_ID = $fields[2]; 
		$flag = $fields[1];
		$map_pos = $fields[3];
		$cigar = $fields[5];
		# stampy reports the mapping position for the first mapping
		# base in the read. Some reads are mapped with an insertion
		# with respect to the reference at the beginning of the read.
		# So the mapping position on a non-diverged reference sequence
		# would be the here reported mapping position minus the insertion.
		# Without this correction, reads with insertions at the beginning
		# of the read would not be detected as coming from the same SbfI restriction site
		# as another read mapped on the opposite strand.
		if($cigar =~ /^(\d+)I/){
			$map_pos -= $1;
		}
		# if the read maps to the reverse strand
		# store the mapping position under reverse reads
		if($flag & 16){ push @{ $map_pos{R} }, $map_pos; }
		else{ push @{ $map_pos{F} }, $map_pos; }
		# while new SAM records can be read in
		while($new_SAM_record = <$IN>){
			# split the record into fields
			@fields = split('\t', $new_SAM_record);
			$flag = $fields[1];
			$cigar = $fields[5];
			$map_pos = $fields[3];
			# if the SAM record belongs to the next contig
			if($fields[2] ne $contig_ID){
				# The single end reads should be 46 base pairs long. If two SE reads come form the 
				# same SbfI restriction site and one read maps on the reverse strand, then the mapping
				# position of the read mapping on the forward strand should be equal to the mapping
				# position of the reverse mapping read + 45 - 3. If that condition is met, store the contig name. 
				# Otherwise, nothing is stored and the next contig is examined.
				MAPPOS: 
				foreach my $pos1 (@{ $map_pos{F} }){
					foreach my $pos2 (@{ $map_pos{R} }){
						if(($pos2 + 45 - 3) == $pos1){ 
							push @linked_RADtag_contigs, $contig_ID;
							last MAPPOS;
						}
					}
				}
				next CONTIG;
			}else{
				if($cigar =~ /^(\d+)I/){
					$map_pos -= $1;
				}
				# store mapping positions
				if($flag & 16){ push @{ $map_pos{R} }, $map_pos; }
				else{ push @{ $map_pos{F} }, $map_pos; }
			}
		}# -- END of File
		# the end of the input file has been reached, thus also the end
		# of records for the last contig. So one last time, check for linked
		# RADtags
		MAPPOS: 
		foreach my $pos1 (@{ $map_pos{F} }){
			foreach my $pos2 (@{ $map_pos{R} }){
				if(abs($pos1 - $pos2) == 41){ 
					push @linked_RADtag_contigs, $contig_ID;
					last MAPPOS;
				}
			}
		}
	}# -- END of CONTIG


	foreach my $contig (@linked_RADtag_contigs){
		print $contig, "\n";
	}

	close $IN;
}


# 
# 
# sub new_SAM_record{
# 
# 	my ($fh) = @_;
# 	my $line = <$fh>;
# 	if(!$line){ return; }
# 	else{ 
# 		my @fields = split('\s+', $line);
# 		return \@fields;
# 	}
# 
# }
