#!/usr/bin/python

# run QC on the samples specified in a given config file
# this is intended to be run from the top level directory for those samples
# we assume that the config file given on the command line points properly to the ipod, chip, and inp directories that are contained under here
# in addition, each of those directories should have a properly filled out read_manifest.txt, and must have raw/, aligned/, and bootstrap/ directories 
# already set up

import sys
import subprocess
import toml
import os
import os.path

SRC_DIR=os.path.dirname( os.path.realpath( __file__ ) )
QC_CMD = "python %s read_manifest.txt" % (os.path.join(SRC_DIR, "run_fastqc_onedir.py"))
CHIP_CMD = "python %s" % (os.path.join(SRC_DIR, "run_chipqc_onedir.py"))
BASE_DIR = os.getcwd()

conf_file = sys.argv[1]
conf_dict = toml.load(conf_file)

ipod_dir = conf_dict["ipod"]["directory"]
chip_dir = conf_dict["chip"]["directory"]
inp_dir = conf_dict["inp"]["directory"]

n_errors = 0

print "Now running fastqc on all preprocessed data..."
print "Working on IPOD samples..."
try:
    os.chdir(ipod_dir)
    subprocess.check_call(QC_CMD, shell=True)
    subprocess.check_call(CHIP_CMD, shell=True)
    os.chdir(BASE_DIR)
except:
    n_errors += 1
    print "Warning: Encountered an error processing IPOD samples"

print "Working on ChIP samples..."
try:
    os.chdir(chip_dir)
    subprocess.check_call(QC_CMD, shell=True)
    subprocess.check_call(CHIP_CMD, shell=True)
    os.chdir(BASE_DIR)
except:
    n_errors += 1
    print "Warning: Encountered an error processing ChIP samples"

print "Working on input samples..."
try:
    os.chdir(inp_dir)
    subprocess.check_call(QC_CMD, shell=True)
    subprocess.check_call(CHIP_CMD, shell=True)
    os.chdir(BASE_DIR)
except:
    n_errors += 1
    print "Warning: Encountered an error processing input samples"

print "Finished with all QC. Encountered %i errors" % n_errors
