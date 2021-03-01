#!/usr/bin/python

# script to run fastqc on each input file in my sequencing manifest
# these should be inspected to make sure everything looks ok, likely before doing ANYTHING else
# by default we run this on the clipped/trimmed reads that are ready to be aligned

# we require fastqc to be in the current path; this was tested with v0.10.1

import subprocess
import shlex
import tempfile
import os
from multiprocessing import Pool
import argparse

# useful constants
READDIR="aligned" ;# directory containing the trimmed reads and fastqc output

F_SUFFIX="_trim_fwd_paired.gz"
R_SUFFIX="_trim_rev_paired.gz"


def run_fastqc_fqfiles(infile_f, infile_r):
  cmdline = "fastqc %s %s --noextract -t 12"  % (infile_f, infile_r)
  subprocess.call(cmdline,shell=True)

intab = open('read_manifest.txt')
samp_prefixes = []
for line in intab:
  this_prefix = line.rstrip().split()[4]
  samp_prefixes.append(this_prefix)

for inprefix in samp_prefixes:
  infile_fwd = os.path.join(READDIR,inprefix + F_SUFFIX)
  infile_rev = os.path.join(READDIR,inprefix + R_SUFFIX)
#  try:
#    test1 = os.stat(infile_fwd)
#    test1 = os.stat(infile_rev) 
#  except:
#    print inprefix
#    print infile_fwd
#    print infile_rev
#    print 'found an error'
#  print '---'
  run_fastqc_fqfiles(infile_fwd,infile_rev)
