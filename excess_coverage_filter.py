#!/usr/bin/env python


# ----------------------------------------
# read in command line arguments
# ----------------------------------------
import argparse

parser = argparse.ArgumentParser(description="Takes BED file, a percentile and list of BAM files and produces a new BED file that has contigs with excess SE read coverage removed. \
        Only SE reads with mapping position 2 are considered. \
        The coverage distributions do not include the count category 0. \
        Filtering is done for individual and across sample coverage for the specified coverage percentile. \
        Contigs with no SE read mapping at position 2 are also automatically removed. \
        Prints the filtered BED file to standard output.")

parser.add_argument("-b", "--bed", help="bed file name including path", type=str)
parser.add_argument("-p", "--percentile", help="percentile of the coverage dist.'s [default: 99.0]", type=float, default=99.0)
parser.add_argument("bamfiles", metavar='BAM files', type=str, nargs='+', help="list of BAM file names, one for each individual")

args = parser.parse_args()


# ----------------------------------------
# load requires modules
# ----------------------------------------
import subprocess32 as sp
import numpy as np
from collections import defaultdict


# ----------------------------------------
# initialise data structures
# ----------------------------------------
# set of contig ids's from contigs with excess individual coverage:
excess_ind_cov = set()
# stores contig id's with excess across sample coverage:
excess_global_cov = set()
# stores individual coverage:
ind_sample_cov = dict([(file, defaultdict(int)) for file in args.bamfiles])
# stores across sample coverage:
across_sample_cov = defaultdict(int)
# stores 99th coverage percentile for each individual:
Q_ind = {}
# list of contigs with SE read mappings from any individual:
all_contigs = set() # contigs with only PE read mappings are completely ignored


# -----------------------------
# get coverage distributions
# -----------------------------
for file in args.bamfiles:
    cmd = "samtools view -f 64 %s | gawk '$4==2' | cut -f 3 | uniq -c" % (file)
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE).stdout
    for line in pipe:
        cov, contig = line.strip().split()
        ind_sample_cov[file][contig] = int(cov)
        Q_ind[file] = 0 # initialise dict
        all_contigs.add(contig)
        across_sample_cov[contig] += int(cov)
    pipe.close()

# ----------------------------------------
# get percentiles of coverage dist.'s
# ----------------------------------------
# get 99th percentiles of individual coverage distributions
for file in args.bamfiles:
    Q_ind[file] = np.percentile(ind_sample_cov[file].values(), args.percentile, interpolation='lower')

# get 99th percentile of across sample coverage distribution
Q_across_sample = np.percentile(across_sample_cov.values(), args.percentile, interpolation='lower')


# ----------------------------------------------------------
# check contigs for excess global, then individual coverage
# ----------------------------------------------------------
for contig in all_contigs:
    if across_sample_cov[contig] > Q_across_sample:
        excess_global_cov.add(contig)
        continue
    for file in args.bamfiles:
        if ind_sample_cov[file][contig] > Q_ind[file]:
            excess_ind_cov.add(contig)
            break


# -------------------------------------------------------
# get new set of contigs that passed filters
# -------------------------------------------------------
# Return a new set with elements in the set 'all_contigs' 
# that are not in the sets 'excess_global_cov' or 'excess_ind_cov':
keep_contigs = all_contigs.difference(excess_global_cov, excess_ind_cov)

bed_in = open(args.bed)

for line in bed_in:
    contig, _, _ = line.strip().split()
    if contig in keep_contigs:
        print line
        
bed_in.close()



