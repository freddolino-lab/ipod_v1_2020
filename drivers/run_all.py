#!/usr/bin/python

# actually run all of the needed commands in this directory

import subprocess
import sys
import os 
import os.path

BASEDIR=os.getcwd()
SRCDIR=os.path.join(os.path.dirname( os.path.realpath( __file__ ) ), "..", "src")

# define the commands that we use for each step

## the following command runs preprocessing and alignment
## it requires one argument: the config file in the working directory with detailed sample information
PR_AL_CMD = "python %s/run_all_alignments.py %%s" % SRCDIR
BS_CMD = "python %s/run_all_bootstraps.py %%s" % SRCDIR
QUANT_CMD = "python %s/analyze_bootstraps/process_bs_files.py %%s" % SRCDIR
QC_CMD = "python %s/run_all_qc.py %%s" % SRCDIR
QNORM_CMD = "python %s/analyze_bootstraps/qnorm_bs_files.py %%s" % SRCDIR

# first just read the set of files to be acted upon
#TARGETS_FILE="all_conditions.txt"
TARGETS_FILE=sys.argv[1]

instr=open(TARGETS_FILE)

all_dirs = []
all_confs = []

for line in instr:
    linearr=line.rstrip().split()
    all_dirs.append(linearr[0])
    all_confs.append(linearr[1])

# now do the first step - preprocessing and alignment

print "Beginning preprocessing and alignment stage..."
print "=============================================="

for dirname,confname in zip(all_dirs, all_confs):
    print ''
    print 'Working on sample name %s' % dirname
    print '----------------------------------'

    os.chdir(dirname)
    print PR_AL_CMD % confname
    subprocess.call(PR_AL_CMD % confname, shell=True)

    print 'done processing sample'
    print '____________________________________'

    os.chdir(BASEDIR)



print "Beginning bootstrapping and QC stage..."
print "=============================================="


for dirname,confname in zip(all_dirs, all_confs):
    print ''
    print 'Working on sample name %s' % dirname
    print '----------------------------------'

    # run the quantitation commands here, making any missing directories of needed
    os.chdir(dirname)
    subprocess.call(BS_CMD % confname, shell=True)
    subprocess.call(QC_CMD % confname, shell=True)

    print '____________________________________'

    os.chdir(BASEDIR)

# quantile normalize all of the samples
print "Doing quantile normalization"
print "=============================================="
for dirname,confname in zip(all_dirs, all_confs):
    print ''
    print 'Working on sample name %s' % dirname
    print '----------------------------------'
    #print QNORM_CMD % conf_files_full
    os.chdir(dirname)
    subprocess.call(QNORM_CMD % confname, shell=True)
    os.chdir(BASEDIR)

print "\nEntering quantitation step"
print "=============================================="

for dirname,confname in zip(all_dirs, all_confs):
    print ''
    print 'Working on sample name %s' % dirname
    print '----------------------------------'

    # run the quantitation commands here, making any missing directories of needed
    os.chdir(dirname)
    print "== running quant cmd"
    print >> sys.stderr, "== running quant cmd"
    print QUANT_CMD % confname
    subprocess.call(QUANT_CMD % confname, shell=True)

    print 'done processing sample'
    print '____________________________________'

    os.chdir(BASEDIR)

print "=============================================="
print "FINISHED WITH ALL SAMPLES"
