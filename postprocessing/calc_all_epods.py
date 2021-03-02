#!/usr/bin/python

# calculate all EPODs for a set of IPOD samples that I am interested in
# we read all of the instances to look at from a table of name/config file pairs

from epod_accessories import *
import sys
import argparse
import os
import os.path
import toml
import subprocess
from multiprocessing import Pool

# set up to parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument('--basedir',help="root directory to find all files of interest (default: current directory)", default=os.getcwd())
parser.add_argument('--outdir',help="root directory for output files (default: current directory)", default=os.getcwd())
parser.add_argument('--numproc',help="number of processors to use (default: 12)", default=12,type=int)
parser.add_argument('--debug',help="run in debug mode (disables parallel processing)", default=False,action='store_true')
parser.add_argument('cond_list', help='File containing a list of condition directories and config files', metavar='condition_list')

args=parser.parse_args()

# define a helper function to actually call epods for a single condition

def do_epod_calls(gr_file_in, outprefix):
    """
    Do all epod calling for a given gr file, writing the results to out_prefix
    """

    print "started EPOD processing for %s" % outprefix

    med512_file = "%s_median512.gr" % outprefix
    med256_file = "%s_median256.gr" % outprefix
    output_epod_file = "%s_epods_v3.gr" % outprefix
    output_peak_file = "%s_epod_locs_v3.txt" % outprefix
    output_epod_file_strict = "%s_epods_strict_v3.gr" % outprefix
    output_peak_file_strict = "%s_epod_locs_strict_v3.txt" % outprefix

    do_runningavg_opt(gr_file_in, med512_file, width=513,genomelength=4641652)
    do_runningavg_opt(gr_file_in, med256_file, width=257,genomelength=4641652)
    identify_epods_v3(med512_file, med256_file, 1024, output_epod_file,delta=25)
    identify_epods_v3(med512_file, med256_file, 1024, output_epod_file_strict,delta=10)

    print "finished EPOD processing for %s" % outprefix


## now actually go through the conditions of interest and run the analysis

if __name__ == "__main__":

    myprocs=Pool(args.numproc)

    conf_str=open(args.cond_list)
    for line in conf_str:
        dirname,conffile = line.rstrip().split()
        cond_conf = toml.load(os.path.join(args.basedir, dirname, conffile))
        gr_file_in = os.path.join(args.basedir, dirname, cond_conf["general"]["output_path"], cond_conf["general"]["out_prefix"] + "_v6rz_chipsub.gr")
        outprefix = os.path.join( args.outdir, cond_conf["general"]["out_prefix"] + "_v6rz_chipsub_epods")

        if args.debug:
            apply(do_epod_calls, [gr_file_in, outprefix])
        else:
            myprocs.apply_async(do_epod_calls, [gr_file_in, outprefix])

    conf_str.close()
    myprocs.close()
    myprocs.join()
