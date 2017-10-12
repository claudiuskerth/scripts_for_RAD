#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="folds 1D spectrum, see eq. 1.2 in Wakeley2009, reads spectrum from stdin, \
        assumes that SFS starts with count for invariant sites")
 
args = parser.parse_args()

from sys import stdin,stdout,stderr

# print >> stdout, stdin.readline().rstrip()

# read in 1D SFS from STDIN, assumes SFS is on one line
sfs = stdin.readline().rstrip().split()

# print '\t'.join(sfs)

n = len(sfs)-1 # sample size in number of alleles per locus

sfs_folded = [0] * (n/2+1)

for i in range(n/2+1):
    sfs_folded[i] = float(sfs[i]) + float(sfs[n-i])
    if i == (n-i):
        sfs_folded[i] = sfs_folded[i]/2.0 
    # print str(i), " + ", str(n-i)

print '\t'.join([str(count) for count in sfs_folded])
