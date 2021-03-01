#!/usr/bin/python

from bs_rep_utils import *
from qnorm_bs_files import *
import sys
import toml
import os.path


# run quantile normalization on all bootstrap output files
# we do the quantile normalization separately for the ipod, inp, and chip signals
# note that all distributions - both the original and bootstrap - are normalized based
#  on the quantiles from the median of the original distributions for that sample type
#  but keep the original and normalized files for the parsed data
# note that in this version, we act on all different conditions simultaneously
# this may enable more consistent normalization across conditions
# I use the gqnorm instead of qnorm output prefix to indicate this
# 
# the single command line argument to this program should be a comma separated list of
#  the directory-level config files for each of the biological conditions of interest
# we assume that each config file resides in the correct root directory for that condition

# these are just helper functions

def do_qnorm_for_samptype(samptype, all_dirs,all_parsed_confs):
    """
    apply quantile normalization across conditions for a given sample type

    We pull together all of the experimental replicates for a given experiment type across many different conditions, quantile normalize them, and then write them back to their parent directories

    samptype: the name of the sample tpe (must match a section in the sample config files)
    all_dirs: ordered list of the top-level directories containing the data
    all_parsed_confs: parsed instances of toplevel config files from the target directories, in the same order as all_parsed_confs

    """

    # first, go through all directories and generate a list of all of the files that should be pulled together here
    all_orig_files = []
    all_bs_files = []

    for dir_name,parsed_conf in zip(all_dirs,all_parsed_confs):
        subdir_name = parsed_conf[samptype]['directory']
        samp_prefixes = parsed_conf[samptype]['sample_names']
        orig_suffix = parsed_conf['general']['orig_suffix']
        bs_suffix = parsed_conf['general']['bs_suffix']
        bs_dir = parsed_conf['general']['bs_dirname']
        all_orig_files += [os.path.join( dir_name, subdir_name, bs_dir, x + orig_suffix) for x in samp_prefixes]
        all_bs_files += [os.path.join( dir_name, subdir_name, bs_dir, x + bs_suffix) for x in samp_prefixes]

    # now we actually apply the normalization,
    #  usin code directly from qnorm_bs_files.py

    q_norm_files( all_orig_files, all_bs_files, out_tag = '_gqnorm')

# here is where the main program starts

if __name__ == '__main__':

    all_files = sys.argv[1].split(",")

    # load information on all of the config files of interest
    all_parsed_confs = []
    all_dirs = []

    for f in all_files:
        this_dir,this_conf = os.path.split(f)
        all_dirs.append(this_dir)
        this_parsed = toml.load(f)
        all_parsed_confs.append(this_parsed)

    for samptype in ["ipod","inp","chip"]:
        do_qnorm_for_samptype(samptype,all_dirs, all_parsed_confs)




