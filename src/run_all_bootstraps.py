#!/usr/bin/python

# run all bootstrapping needed for ipod, chip, and input samples specified in a given config file
# this is intended to be run from the top level directory for those samples
# we assume that the config file given on the command line points properly to the ipod, chip, and inp directories that are contained under here
# in addition, each of those directories should have a properly filled out read_manifest.txt, and must have raw/, aligned/, and bootstrap/ directories 
# already set up
# we then do a standard set of bootstrapping where we generate both occupancy traces for the original data, and 200 boostrap replicates, at 5 bp resolution

# the following external programs need to be installed and in the PATH:
#  samtools (tested with 0.1.19-96b5f2294a -- newer versions may not work!)

import sys
import subprocess
import toml
import os
import multiprocessing
import tempfile


conf_file = sys.argv[1]
conf_dict = toml.load(conf_file)

ipod_dir = conf_dict["ipod"]["directory"]
chip_dir = conf_dict["chip"]["directory"]
inp_dir = conf_dict["inp"]["directory"]

BS_SAM_THREADS = conf_dict["bootstrap"]["samtools_threads"]
BS_BOOT_THREADS = conf_dict["bootstrap"]["bootstrap_threads"]

# path information needed for bootstrapping program
SRCDIR = os.path.dirname( os.path.realpath( __file__ ) )
BINDIR=os.path.join( SRCDIR, "bootstrapping")
PARSE_CMD = "python %s/bootstrap_sam_file.py parse --paired %s %s"
SAMPLE_CMD = "python %s/bootstrap_sam_file.py  sample %s.npy %s.npy 4641652 --num_samples 200 --resolution 5"
ORIG_CMD = "python %s/bootstrap_sam_file.py  sample %s.npy %s.npy 4641652 --identity --resolution 5"


n_errors = 0

# here we define a helper function to do the bootstrap part ONLY

def preprocess_bootstrap( bamfile, outprefix, nthreads ):
    # do the preprocessing necessary to set up a bam file for bootstrapping
    # this follows the procedure that I had established in /data/petefred/st_lab_work/ipod_hr/bin/bootstrapping/preprocess_bamfile.py
    # this does the initial setup and parse stages, but not the actual bootstrap resampling

    # first get a sorted sam file as needed by the parse command
    tmp_file = tempfile.NamedTemporaryFile()
    filter_cmd = "samtools sort -n -o -@ %i %s %s_tmp | samtools view -f 3 -F 2308 -q 30 - | sed s/\"#0\/4\t\"/\"#0\t\"/g> %s" % (nthreads, bamfile, bamfile, tmp_file.name)
    print filter_cmd
    subprocess.call(filter_cmd, shell=True)

    # now run the parse step
    print (PARSE_CMD % (BINDIR,tmp_file.name, outprefix+"_parsed"))
    subprocess.call(PARSE_CMD % (BINDIR,tmp_file.name, outprefix+"_parsed"),shell=True)

    # and clean up
    tmp_file.close()


def do_bootstrap( outprefix ):
    # assuming that parsing has already completed for a specified bam file, now we actually run the needed resampling

    # do the actual bootstrapping 
    print (SAMPLE_CMD % (BINDIR,outprefix+"_parsed", outprefix+"_resampled"))
    subprocess.call(SAMPLE_CMD % (BINDIR,outprefix+"_parsed", outprefix+"_resampled"),shell=True)

    # also generate a file with the original occupancies
    print (ORIG_CMD % (BINDIR,outprefix+"_parsed", outprefix+"_original"))
    subprocess.call(ORIG_CMD % (BINDIR,outprefix+"_parsed", outprefix+"_original"),shell=True)

# now set up the pool that I will use for most of the work
bs_pool = multiprocessing.Pool(BS_BOOT_THREADS)
all_thr = [] ;# this is a list that will contain all of the result objects

print "Beginning work on bootstrapping"
print "Working on IPOD samples..."

ipod_prefixes = conf_dict["ipod"]["sample_names"]
bam_files = [ os.path.join( ipod_dir, "aligned", pref + "_bowtie2_sorted.bam" ) for pref in ipod_prefixes ]
out_names = [ os.path.join( ipod_dir, "bootstrap", pref ) for pref in ipod_prefixes ]

# make the bootstrap directory if it does not already exist
if not(os.path.isdir(os.path.join( ipod_dir, "bootstrap"))):
    os.mkdir(os.path.join( ipod_dir, "bootstrap"))

for bamname, outpref in zip(bam_files, out_names):
    print "Preprocessing %s" % bamname

    try:
        preprocess_bootstrap( bamname, outpref, BS_SAM_THREADS)
        print "Background processing bootstrap for %s" % bamname
        all_thr.append(bs_pool.apply_async( do_bootstrap, [outpref] ))
        #subprocess.check_call(run_cmd, shell=True)
    except:
        n_errors += 1
        print "Warning: Encountered an error processing IPOD sample %s" % bamname


print "Working on ChIP samples..."

chip_prefixes = conf_dict["chip"]["sample_names"]
bam_files = [ os.path.join( chip_dir, "aligned", pref + "_bowtie2_sorted.bam" ) for pref in chip_prefixes ]
out_names = [ os.path.join( chip_dir, "bootstrap", pref ) for pref in chip_prefixes ]

# make the bootstrap directory if it does not already exist
if not(os.path.isdir(os.path.join( chip_dir, "bootstrap"))):
    os.mkdir(os.path.join( chip_dir, "bootstrap"))

for bamname, outpref in zip(bam_files, out_names):
    print "Preprocessing %s" % bamname

    try:
        preprocess_bootstrap( bamname, outpref, BS_SAM_THREADS)
        print "Background processing bootstrap for %s" % bamname
        all_thr.append(bs_pool.apply_async( do_bootstrap, [outpref] ))
    except:
        n_errors += 1
        print "Warning: Encountered an error processing ChIP sample %s" % bamname


print "Working on Input samples..."

inp_prefixes = conf_dict["inp"]["sample_names"]
bam_files = [ os.path.join( inp_dir, "aligned", pref + "_bowtie2_sorted.bam" ) for pref in inp_prefixes ]
out_names = [ os.path.join( inp_dir, "bootstrap", pref ) for pref in inp_prefixes ]

# make the bootstrap directory if it does not already exist
if not(os.path.isdir(os.path.join( inp_dir, "bootstrap"))):
    os.mkdir(os.path.join( inp_dir, "bootstrap"))

for bamname, outpref in zip(bam_files, out_names):
    print "Preprocessing %s" % bamname

    try:
        preprocess_bootstrap( bamname, outpref, BS_SAM_THREADS)
        print "Background processing bootstrap for %s" % bamname
        all_thr.append(bs_pool.apply_async( do_bootstrap, [outpref] ))
    except:
        n_errors += 1
        print "Warning: Encountered an error processing Input sample %s" % bamname


# now collect all of the threads
bs_pool.close()
bs_pool.join()

n_pool_err = 0
for res in all_thr:
    if not res.successful():
        n_pool_err += 1

print "Finished with bootstrapping. Encountered %i errors with preprocessing and %i errors during the bootstrap portion" % (n_errors, n_pool_err)
