#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: cluster.pl
#
#        USAGE: ./cluster.pl  
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
#      CREATED: 08/03/15 12:10:40
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

# TODO: group reads by SbfI position, then cluster within groups based on sequence upstream of SbfI

my (%cluster, $pos, %pos_group, $mismatch_count, $merged, $pre_SbfI_seq_1, $pre_SbfI_seq_2, $threshold);

while(<>){
	chomp;
	if(/CCTGCAGG/){ 
		$pos = $-[0]+1;
		push @{ $pos_group{$pos} }, $_;
	}
}

foreach my $pos (sort {$a <=> $b} keys %pos_group){
	$threshold = ($pos - 1 - 6)/6;
	my %read_group_cluster;
	foreach my $read ( @{ $pos_group{$pos} } ){
		($pre_SbfI_seq_1) = $read =~ m/^TGCAGG(.*)CCTGCAGG.*/;
		$merged = 0;
		CLUSTER: foreach my $cluster_rep (keys %read_group_cluster){
			($pre_SbfI_seq_2) = $cluster_rep =~ m/^TGCAGG(.*)CCTGCAGG.*/;
			$mismatch_count = ($pre_SbfI_seq_2 ^ $pre_SbfI_seq_1) =~ tr/\001-\255//;
			if($mismatch_count <= $threshold){
				push @{ $read_group_cluster{$cluster_rep} }, $read;
				$merged = 1;
				last CLUSTER;
			}
			foreach my $seq ( @{ $read_group_cluster{$cluster_rep} } ){
				($pre_SbfI_seq_2) = $seq =~ m/^TGCAGG(.*)CCTGCAGG.*/;
				$mismatch_count = ($pre_SbfI_seq_2 ^ $pre_SbfI_seq_1) =~ tr/\001-\255//;
				if($mismatch_count <= $threshold){
					push @{ $read_group_cluster{$cluster_rep} }, $read;
					$merged = 1;
					last CLUSTER;
				}
			}
		}
		unless($merged){
			$read_group_cluster{$read} = [];
		}
	}
	# a bit more transitive cluster merging
	# this is not perfect but probably good enough at the moment 
	my @cluster_reps = sort keys %read_group_cluster;
	REP: for(my $i = 0; $i < @cluster_reps; $i++){
		($pre_SbfI_seq_1) = $cluster_reps[$i] =~ m/^TGCAGG(.*)CCTGCAGG.*/;
		for(my $o = $i+1; $o < @cluster_reps; $o++){
			foreach my $seq ( @{ $read_group_cluster{$cluster_reps[$o]} } ){
				($pre_SbfI_seq_2) = $seq =~ m/^TGCAGG(.*)CCTGCAGG.*/;
				$mismatch_count = ($pre_SbfI_seq_2 ^ $pre_SbfI_seq_1) =~ tr/\001-\255//;
				if($mismatch_count <= $threshold){
					push @{ $read_group_cluster{$cluster_reps[$o]} }, $cluster_reps[$i];
					push @{ $read_group_cluster{$cluster_reps[$o]} }, @{ $read_group_cluster{$cluster_reps[$i]} };
					delete $read_group_cluster{$cluster_reps[$i]};
					next REP;
				}
			}
		}
	}
	$cluster{$pos} = \%read_group_cluster;
}

foreach my $pos (sort {$a <=> $b} keys %cluster){
	foreach my $cluster_rep (sort keys %{ $cluster{$pos} }){
		print $cluster_rep, "\n";
		if(scalar @{ $cluster{$pos}->{$cluster_rep} }){
			print join("\n", @{ $cluster{$pos}->{$cluster_rep} }), "\n";
		}
		print "\n";
	}
}
