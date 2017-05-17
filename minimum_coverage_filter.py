#!/usr/bin/env python

# ----------------------------------------
# read in command line arguments
# ----------------------------------------
import argparse

parser = argparse.ArgumentParser(description="Reads depth table from samtools depth command via STDIN and prints to STDOUT sites with \
        minimum coverage in a minimum number of individuals.")

#parser.add_argument("-d", "--depthfile", help="file name of output from samtools depth command", type=str)
parser.add_argument("-mc", "--minimum_coverage", help="minimum fold coverage that is required in at least --minimum_individual individuals (default 1)", type=int, default=1)
parser.add_argument("-mi", "--minimum_individual", help="minimum number of individuals with at least --minimum_coverage coverage (default 15)", type=int, default=15)

args = parser.parse_args()

# print args.minimum_coverage
# print args.minimum_individual

# import gzip
# 
# depth_fh = gzip.open(args.depthfile)

import sys

import numpy as np

for line in sys.stdin:
    fields = line.strip().split()
    depths = np.array(fields[2:], dtype='int8')
    if sum( depths >= args.minimum_coverage ) >= args.minimum_individual:
        print line # print to standard output the contig name, position and depths

# depth_fh.close()
