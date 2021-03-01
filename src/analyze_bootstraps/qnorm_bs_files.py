#!/usr/bin/python

import numpy
import scipy.stats
from bs_rep_utils import *
import toml
import os.path
import sys


# run quantile normalization on all bootstrap output files
# we do the quantile normalization separately for the ipod, inp, and chip signals
# note that all distributions - both the original and bootstrap - are normalized based
#  on the quantiles from the median of the original distributions for that sample type
# we delete the un-normalized bootstrap replicate matrices (adding a suffix to the name),
#  but keep the original and normalized files for the parsed data


def calc_qnorm_base( input_arrays ):
    # calculate the averaged quantile-wise distribution over a bunch of targets
    # we take the median of the input values at each quantile to define that quantile
    # note that we assume that there are no missing values, and that all inputs
    # have the same length


    # first we sort the inputs column-wise
    sorted_arrs = []
    for i in range(len(input_arrays)):
        sorted_arrs.append( numpy.sort(input_arrays[i],kind='heapsort') )

    # then pull those sorted arrays together and figure out the median quantile distribution
    sorted_mat = numpy.vstack( sorted_arrs )
    median_qs = numpy.median(sorted_mat, axis=0)

    return median_qs

def q_norm_vec( input_vals, target_vals ):
    # quantile normalize input_vals to match the distribution in target_vals

    rank_vec = scipy.stats.rankdata(input_vals, method='ordinal')
    return target_vals[rank_vec - 1]

def qnorm_bootstrap_mat(fn, target_dist, outfile=None):
    # apply bootstrap normalization separately to each column in
    #  a matrix, using target_dist as the target set of values
    # if outfile is None, we overwrite the input

    full_mat = numpy.load(fn)

    nvals,nboot=full_mat.shape

    for b in range(nboot):
        this_vec = full_mat[:,b]
        full_mat[:,b] = q_norm_vec(this_vec, target_dist)

    if outfile is None:
        numpy.save(fn,full_mat)
    else:
        numpy.save(outfile, full_mat)

def q_norm_files(orig_files, bs_files, out_tag="_qnorm"):
    # read and quantile normalize a set of files
    # we read all of the members of orig_files, calculate the target distribution based
    #  on them, and then write the normalized version of each
    # Then, every bootstrap replicate in each of the files in bs_files is normalized 
    #  TO THE SAME TARGET DISTRIBUTION
    # In the process we over-write and then delete the bootstrap matrices
    #  we do retain the original data both before and after normalization, though

    # first calculate the target distribution and rewrite the original files
    orig_vecs = [ numpy.load( f ).ravel() for f in orig_files ]

    print 'calculating target distribution.......'
    target_distr = calc_qnorm_base(orig_vecs)

    print 'normalizing original data......'
    qnorm_vecs = [ q_norm_vec( v, target_distr) for v in orig_vecs ]
    for f, q in zip(orig_files, qnorm_vecs):
        o_pref,o_suff = os.path.splitext(f)
        #o_name = '_qnorm' + f
        o_name= o_pref + out_tag + o_suff
        numpy.save( o_name, q)

    # now do the same for each of the bootstrap files
    print 'normalizing bootstrap data.....'
    for fn in bs_files:
        o_pref,o_suff = os.path.splitext(fn)
        o_name = o_pref + out_tag + o_suff
        print "normalizing %s to %s" % (fn,o_name)
        qnorm_bootstrap_mat(fn,target_distr, outfile=o_name)
        #o_name =  '_qnorm' + f
        # use this instead of the shutil above to remove old file os.rename( fn, o_name )

# here is where the main program starts

if __name__ == '__main__':

    conf_file = sys.argv[1]
    conf_dict = toml.load(conf_file)

    # figure out some global parameters
    bs_suffix = conf_dict['general']['bs_suffix']
    bs_dir = conf_dict['general']['bs_dirname']
    orig_suffix = conf_dict['general']['orig_suffix']
    out_prefix = os.path.join(conf_dict['general']['output_path'],conf_dict['general']['out_prefix'])

    # run quantile normalization on each set of samples

    ## first the ipod samples
    ipod_dir = conf_dict['ipod']['directory']
    ipod_prefixes = conf_dict['ipod']['sample_names']
    ipod_orig_files = [ os.path.join( ipod_dir, bs_dir, x + orig_suffix) for x in ipod_prefixes]
    ipod_bs_files = [ os.path.join( ipod_dir, bs_dir, x + bs_suffix) for x in ipod_prefixes]
    q_norm_files( ipod_orig_files, ipod_bs_files )

    ## next the input samples
    inp_dir = conf_dict['inp']['directory']
    inp_prefixes = conf_dict['inp']['sample_names']
    inp_orig_files = [ os.path.join( inp_dir, bs_dir, x + orig_suffix) for x in inp_prefixes]
    inp_bs_files = [ os.path.join( inp_dir, bs_dir, x + bs_suffix) for x in inp_prefixes]
    q_norm_files( inp_orig_files, inp_bs_files )

    ## and finally the chip samples
    chip_dir = conf_dict['chip']['directory']
    chip_prefixes = conf_dict['chip']['sample_names']
    chip_orig_files = [ os.path.join( chip_dir, bs_dir, x + orig_suffix) for x in chip_prefixes]
    chip_bs_files = [ os.path.join( chip_dir, bs_dir, x + bs_suffix) for x in chip_prefixes]
    q_norm_files( chip_orig_files, chip_bs_files )




