#!/usr/bin/env python

from __future__ import division
# ----------------------------
# parse command line arguments
# ----------------------------
import argparse
parser = argparse.ArgumentParser(description="apply Ludovic's correction to a bunch of unfolded spectra in dadi input format, the degree of correction is determined individually for each SFS")
parser.add_argument('sfs_file_names', nargs='+')
args = parser.parse_args()


# ---------------
# import modules
# ---------------
import sys
sys.path.insert(0, "/home/claudius/Downloads/dadi/")
import dadi
import numpy as np
# turn on floating point division by default, old behaviour via '//'
from scipy.optimize import minimize_scalar


# -----------------
# define functions
# -----------------
# define function that modifies a SFS with a given p
def modify(p, sfs):
    """
    modifies SFS with given p
    subtracts proportion p from even count classes and adds to corresponding 
    haploid count class, e. g. from 2 onto 1 or from 16 onto 8
    returns modified sfs
    """
    fs = sfs.copy() # it is necessary to make sure that a copy of the observed SFS is
    # made here for each function evaluation

    # modify observed SFS with given p:
    for i in range(len(fs)):
        if not fs.mask[i]:
            if not i%2:
                deduct = fs[i] * p
                fs[i] -= deduct
                j = i//2
                fs[j] += deduct
    return fs


# define cost function used for optimisation
def f(p, snm, sfs):
    """
    p: proportion of all SNP's on the X chromosome [float, 0<p<1]
    snm: standard neutral model spectrum (optimally scaled)
    sfs: observed SFS 
    """
    # modify sfs
    fs = modify(p, sfs)

    # return sum of squared deviations of modified SFS with snm spectrum:
    return np.sum( (fs - snm)**2 )  # note, due to the mask it is necessary to use the numpy sum function 


def get_snm(fs):
    """
    takes folded SFS and returns optimally scaled standard neutral model specrum
    """
    # make the extrapolating version of the standard neutral model function
    func_ex = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.snm)

    # setting the smallest grid size slightly larger than the largest population sample size
    pts_l = [40, 50, 60]

    ns = fs.sample_sizes

    # calculate unfolded AFS under standard neutral model (up to a scaling factor theta)
    snm = func_ex(0, ns, pts_l)
    snm = snm.fold()

    # scaled snm spectrum for fs
    return dadi.Inference.optimally_scaled_sfs(snm, fs)


# -----------------
# main programme
# -----------------
if __name__ == '__main__':
    sfs = []

    for fname in args.sfs_file_names:
        # print fname
        with open(fname, "r") as fh:
            # get SFS from file
            sfs = dadi.Spectrum.from_file(fh)
            sfs = sfs.fold()
            # get optimally scaled standard neutral model
            snm = get_snm(sfs)
            # optimise p
            res = minimize_scalar(f, method = 'bounded', bounds = (0.0, 1.0), args = (snm, sfs))
            # save modified spectrum to file
            dadi.Spectrum.tofile(modify(res.x, sfs), fid=fname + ".corr", comment_lines=['Ludovics correction applied, p=' + str(round(res.x, 2))])


