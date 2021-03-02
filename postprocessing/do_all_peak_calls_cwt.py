#!/usr/bin/python

# calculate peak calls for all of my ipod samples of interest
# we read all of the instances to look at from a table of name/config file pairs
# here we do the calls using the function find_peaks_cwt from scipy, at a variety of thresholds

import argparse
import os
import os.path
import toml
import subprocess
import multiprocessing


# the following command takes three arguments: an input .gr file, an output .gr file, and a threshold value for peak calls
SRC_DIR = os.path.dirname( os.path.realpath( __file__ ) )
PEAK_CALL_SCRIPT="python %s/../src/call_peaks_cwt.py %%s %%s %%f 6" % SRC_DIR

# set up to parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument('--basedir',help="root directory to find all files of interest (default: current directory)", default=os.getcwd())
parser.add_argument('--outdir',help="root directory for output files (default: current directory)", default=os.getcwd())
parser.add_argument('--scoretype',help="type of score files to use for peak calling (default: _v6rzlog10p_chipsub )", default="_v6rzlog10p_chipsub")
parser.add_argument('--threads', help='Number of threads to use', type=int, default=8)
parser.add_argument('cond_list', help='File containing a list of condition directories and config files', metavar='condition_list')

args=parser.parse_args()

# helper function to be called for a single case
def do_calls_overlaps( infile, outprefix, cutoff):
    print "Working on %s / %s / %f" % (infile, outprefix, cutoff)
    run_cmd = PEAK_CALL_SCRIPT % ( infile, outprefix + "_cutoff%f_cwt_peaks" % cutoff, cutoff)
    subprocess.call(run_cmd, shell=True)

# now go through the conditions of interest and run the analysis
# we actually call the peaks, and then compare them to tfbs lists

conf_str=open(args.cond_list)
job_pool = multiprocessing.Pool(args.threads)


for line in conf_str:
    print line
    dirname,conffile = line.rstrip().split()
    cond_conf = toml.load(os.path.join(args.basedir, dirname, conffile))
    for cutoff in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,7.0,8.0,10.0, 12.0, 0.0, 0.25, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]:
        job_pool.apply_async( do_calls_overlaps, [os.path.join(args.basedir, dirname, cond_conf["general"]["output_path"], cond_conf["general"]["out_prefix"] + args.scoretype + ".gr"), os.path.join(args.outdir, cond_conf["general"]["out_prefix"] + args.scoretype), cutoff])

job_pool.close()
job_pool.join()

