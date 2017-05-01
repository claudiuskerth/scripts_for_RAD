#!/usr/bin/env python

# see the IPython notebook coverage_filter.ipynb for the testing of this code

# ----------------------------------------
# read in command line arguments
# ----------------------------------------
import argparse

parser = argparse.ArgumentParser(description="takes output from samtools depth command and filters sites against excessive and for enough coverage; prints filtered sites to STDOUT")

parser.add_argument("--global_coverage_percentile", help="percentile of the global coverage distribution that the site's global coverage must not exceed (default 99)", type=float, default=99.0)
parser.add_argument("--individual_coverage_percentile", help="percentile of the individual coverage distributions that the site's individual coverages must not exceed in any individual (default 99)", type=float, default=99.0)
parser.add_argument("--minimum_coverage", help="minimum fold coverage that is required in at least --minimum_individual individuals (default 1)", type=int, default=1)
parser.add_argument("--minimum_individual", help="minimum number of individuals with at least --minimum_coverage coverage (default 15)", type=int, default=15)
parser.add_argument("depth_file", help="file name for depth table from samtools depth command", type=str)

args = parser.parse_args()


# ----------------------------
# define percentile function
# ---------------------------
def compute_percentile_from_dict(my_dict, percentile):
    """
    my_dict : dictionary, containing depths as keys and their counts as values
    percentile : float in range of [0,100], percentile to compute, which must be between 0 and 100 inclusive.

    returns : key in my_dict that is less than or equal to the specified percentile of the distribution of "keys"
    (i. e. values are counts) stored in my_dict
    """
    assert percentile < 100 and percentile > 0, "percentile must be a number between 0 and 100"

    # get sum of counts
    total_cnt = sum(my_dict.values())

    perc = float() # initialise percentile of the coverage distribution
    cum_cnt = float()

    for d in sorted(my_dict.keys(), reverse=True): # search from the upper tail of the distribution
        cum_cnt += my_dict[d]
        perc = cum_cnt / total_cnt * 100
        if perc < (100 - percentile):
            continue
        else:
            return d

# ----------------------------
# open depth file
# ---------------------------
import gzip 

depth_fh = gzip.open(args.depth_file)


# ------------------------------
# determine coverage thresholds
# ------------------------------

from collections import defaultdict

# initialisze global coverage count dictionary
global_cov_dist = defaultdict(int)

# initialise individual coverage count dictionaries
fields = depth_fh.readline().rstrip().split('\t')
num_inds = len(fields[2:])
ind_cov_dists = dict([(ind, defaultdict(int)) for ind in range(num_inds)])
depth_fh.seek(0) # seek back to beginning of file

# first read through file
for line in depth_fh:
    fields = line.rstrip().split('\t')
    depths = fields[2:]
    # get global coverage
    across_sample_depth = sum([int(d) for d in depths])
    # count global coverages in dictionary
    global_cov_dist[across_sample_depth] += 1
    # count individual coverages in dictionaries
    for ind, d in enumerate(depths):
        ind_cov_dists[ind][int(d)] += 1 


# set global coverage threshold
global_coverage_threshold = compute_percentile_from_dict(global_cov_dist, args.global_coverage_percentile)
# set individual overage thresholds
ind_coverage_thresholds = [compute_percentile_from_dict(ind_cov_dists[i], args.individual_coverage_percentile) for i in range(num_inds)]


# ------------------------------------------------------------
# filter sites based on coverage thresholds
# ------------------------------------------------------------
#
# second read through file
#
depth_fh.seek(0) # seek back to beginning of file

for line in depth_fh:
    fields = line.rstrip().split('\t')
    depths = fields[2:]
    across_sample_depth = sum([int(d) for d in depths])
    if across_sample_depth > global_coverage_threshold:
        continue
    excess_depth = False
    sufficient_depth = 0
    for ind, d in enumerate(depths):
        if int(d) > ind_coverage_thresholds[ind]:
            excess_depth = True
            break
        elif int(d) >= args.minimum_coverage:
            sufficient_depth += 1
    if not excess_depth and sufficient_depth >= args.minimum_individual:
        print '\t'.join(fields) # prints to STDOUT

depth_fh.close()






















