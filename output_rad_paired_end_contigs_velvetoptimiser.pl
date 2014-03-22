#!/usr/bin/env perl

# output_rad_paired_end_contigs
#
# Take a set of RAD tags and a set of files containing the reads for all
# individuals, get the paired ends for each tag, assemble using phrap and
# output best contig
# Reads in format
# Singleseq Singlequal Pairedseq Pairedqual
# Tags in format
# Tag Count1 Count2 Count3...
# with header specifying individuals
# Author: John Davey john.davey@ed.ac.uk
# Begun 31/3/10
# v0.1: Runs on single contig files
# v0.2:21/06/10  Modified for array jobs
#############################################################################
###                                                                       ###
###                                 CODE                                  ###
###                                                                       ###
#############################################################################

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Data::Dumper;

use Bio::SeqIO;

# Autoflush output so reporting on progress works
$| = 1;

my @pair_filenames = glob('*');

foreach my $tag_header (sort @pair_filenames) {
    my ($tag,$tag_pair_num) = split /_/, $tag_header;

    # Run assembler and get longest contig from output
    #VelvetOptimiser

    my $tag_dir = "$tag_header\_dir";

    mkdir $tag_dir;
    chdir $tag_dir;
    my $read_length = length($tag);

    system(
    "VelvetOptimiser.pl -f \'-fasta -short ../$tag_header\' -o \'-read_trkg yes -amos_file yes -min_contig_lgth $read_length\' -c max  1>/dev/null 2>/dev/null"
    );

    my @contig_files = glob('auto_data_*/contigs.fa');
    chdir("..");

    # If only one contig file exists, as it should if VelvetOptimiser runs
    # successfully, write out the assembled contig(s)
    if ( @contig_files == 1 ) {

        open my $contigs_out_file, ">", "$tag_header\.out";
        print $contigs_out_file ">$tag_header\_tag\n$tag\n";
        my $contigs = new Bio::SeqIO(
            -format => 'FASTA',
            -file   => "$tag_dir/$contig_files[0]"
        );
        my $seq_count = 0;
        while ( my $seq = $contigs->next_seq ) {
            $seq_count++;
            my $seq_out = $seq->seq();
            my $title = $seq->display_id();
            my $length = "";
            my $coverage = "";
            if ($title =~ /^(.+)length_(.+)_cov_(.+)$/) {
                $length = $2;
                $coverage = $3;
            }
            print $contigs_out_file ">$tag_header\_pair_$seq_count\_length\_$length\_cov_$coverage\n$seq_out\n";
        }
        close $contigs_out_file;
    }
    else 
    {
        system("rm -r $tag_dir");
    }
}
