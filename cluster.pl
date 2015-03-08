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


my (%cluster, $mismatch_count, $merged);

while(<>){
	chomp;
	$mismatch_count = 0;
	$merged = 0;
	CLUSTER: foreach my $cluster_rep (keys %cluster){
		$mismatch_count = ($cluster_rep ^ $_) =~ tr/\001-\255//;
		if($mismatch_count <= 5){
			push @{ $cluster{$cluster_rep} }, $_;
			$merged = 1;
			last CLUSTER;
		}
		foreach my $seq ( @{ $cluster{$cluster_rep} } ){
			$mismatch_count = ($seq ^ $_) =~ tr/\001-\255//;
			if($mismatch_count <= 5){
				push @{ $cluster{$cluster_rep} }, $_;
				$merged = 1;
				last CLUSTER;
			}
		}
	}
	unless($merged){
		$cluster{$_} = [];
	}
}

# a bit more transitive cluster merging
# this is not perfect but probably good enough at the moment 
my @cluster_reps = sort keys %cluster;
REP: for(my $i = 0; $i < @cluster_reps; $i++){
	for(my $o = $i+1; $o < @cluster_reps; $o++){
		foreach my $seq ( @{ $cluster{$cluster_reps[$o]} } ){
			$mismatch_count = ($seq ^ $cluster_reps[$i]) =~ tr/\001-\255//;
			if($mismatch_count <= 5){
				push @{ $cluster{$cluster_reps[$o]} }, $cluster_reps[$i];
			   	push @{ $cluster{$cluster_reps[$o]} }, @{ $cluster{$cluster_reps[$i]} };
				delete $cluster{$cluster_reps[$i]};
				next REP;
			}
		}
	}
}

foreach my $cluster_rep (sort keys %cluster){
	print $cluster_rep, "\n";
	if(scalar @{ $cluster{$cluster_rep} }){
		print join("\n", @{ $cluster{$cluster_rep} }), "\n";
	}
	print "\n";
}
