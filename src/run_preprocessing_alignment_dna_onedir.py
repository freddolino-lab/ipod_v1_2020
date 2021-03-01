#!/usr/bin/python

# script for doing my standardized preprocessing and alignment on DNA files
# we do standardized quality trimming, adapter removal, and then alignment to the E. coli MG1655 genome (U00096.3)
# we expect a fairly standardized directory layout:
#  There should be a 'raw' subdirectory containing the raw reads, and an 'aligned' subdirectory to contain the aligned data
#  In addition, the top level directory should have a space-delimited file called 'seq_manifest.txt' 
#  This file needs to have a series of input prefix/output prefix pairs, and will be used to identify which runs to look at

# in addition, the following programs must be present in the current system PATH:
# cutadapt (tested with 1.8.1)
# trimmomatic (tested with 0.32)
# bowtie2 (tested with 2.1.0)

import sys
import subprocess
import shlex
import tempfile
import os
import os.path
import argparse
from multiprocessing import Pool
from bunch import Bunch

NPROC=12
BINDIR = os.path.dirname( os.path.realpath( __file__ ) )
SEQ_DB="%s/../db/u00096_3" % BINDIR
BTBIN="%s/../bin/bowtie2" % BINDIR

## Set up defined constants that should be universal
READDIR='raw'
ALDIR='aligned'

E_COLI_DBFILE = SEQ_DB

# set up the needed directories if they are not already present
if not(os.path.isdir(READDIR)):
    os.mkdir(READDIR)

if not(os.path.isdir(ALDIR)):
    os.mkdir(ALDIR)

# the actual input to this program should be a single space-delmited file
# each line should have, in order:
# Freads Rreads adapseq phredbase outprefix

#now, we go parse that file
all_samples = []

instr=open(sys.argv[1])
for line in instr:
    if line[0] == '#':
        continue

    linearr = line.rstrip().split()
    this_samp = Bunch()
    this_samp.ffile = linearr[0]
    this_samp.rfile = linearr[1]
    this_samp.adapseq = linearr[2]
    this_samp.phredbase = int(linearr[3])
    this_samp.outprefix = linearr[4]
    all_samples.append(this_samp)

# define some functions that will be used in the rest of the script

def preprocess_gz_file(samp):
  # do some initial preprocessing of a gz file, including trimming and quality score filtering
  # we need to give it a bunch that has defined all of the
  # data listed higher up in this file (don't use this function outside of this script!!)

  infile_1 = samp.ffile
  infile_2 = samp.rfile
  outprefix=samp.outprefix
  PHRED_BASE=samp.phredbase
  ADAP_SEQ = samp.adapseq

  if infile_1[-3:] == ".gz":
    DCPROG = 'zcat'
  elif infile_1[-4:] == ".bz2":
    DCPROG = 'bzcat'
  else:
    raise("Could not determine the decompression program to use")

  if DCPROG == "bzcat":
    in1 = tempfile.NamedTemporaryFile(suffix='.fastq')
    in2 = tempfile.NamedTemporaryFile(suffix='.fastq')
    infile_fwd = in1.name
    infile_rev = in2.name
    cmd1="%s %s > %s" % (DCPROG, infile_1, infile_fwd)
    cmd2="%s %s > %s" % (DCPROG, infile_2, infile_rev)
    subprocess.call(cmd1,shell=True)
    subprocess.call(cmd2,shell=True)
 


  else:
    infile_fwd = infile_1
    infile_rev = infile_2

  cutfile_fwd = os.path.join(ALDIR, outprefix+"_fwd_cutadap.fq.gz")
  cutfile_rev = os.path.join(ALDIR, outprefix+"_rev_cutadap.fq.gz")

  #  do some quality trimming and write a processed file


  # first clip the adapter sequences
  cutadapt_cmd = "cutadapt --quality-base=%i -a %s -A %s -n 3 --match-read-wildcards -o %s -p %s %s %s > %s_cutadapt.log 2> %s_cutadapt.err" % ( PHRED_BASE,ADAP_SEQ, ADAP_SEQ, cutfile_fwd, cutfile_rev, infile_fwd, infile_rev, outprefix,outprefix)
  print cutadapt_cmd
  subprocess.call(cutadapt_cmd,shell=True)

  # next do quality trimming -- trim true crap from the 3' end, and then look for a  sliding window of 4 bp with qualities above 15
  # we avoid doing 5' end trimming so we don't move where the true end of the read is
  # drop the read if we have less than 10 bases after this is done
  trim_fwd_paired = os.path.join( ALDIR, outprefix + "_trim_fwd_paired.gz" )
  trim_fwd_unpaired = os.path.join( ALDIR, outprefix + "_trim_fwd_unpaired.gz" )
  trim_rev_paired = os.path.join( ALDIR, outprefix + "_trim_rev_paired.gz" )
  trim_rev_unpaired = os.path.join( ALDIR, outprefix + "_trim_rev_unpaired.gz" )
  trim_cmd = "trimmomatic PE -threads %i -phred%i %s %s %s %s %s %s TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:10 2> %s_trimmomatic.err" % (NPROC, PHRED_BASE,cutfile_fwd, cutfile_rev, trim_fwd_paired, trim_fwd_unpaired, trim_rev_paired, trim_rev_unpaired, outprefix )
  print trim_cmd
  subprocess.call(trim_cmd,shell=True)

  if DCPROG == "bzcat":
    in1.close()
    in2.close()


def run_bowtie(prefix, phredbase, db=E_COLI_DBFILE):
  # run bowtie to align reads to a db
  fwd = os.path.join(ALDIR,prefix+"_trim_fwd_paired.gz")
  rev = os.path.join(ALDIR,prefix+"_trim_rev_paired.gz")
  fwd_unpaired = os.path.join( ALDIR, prefix + "_trim_fwd_unpaired.gz" )
  rev_unpaired = os.path.join( ALDIR, prefix + "_trim_rev_unpaired.gz" )
  samout = os.path.join(ALDIR,prefix+"_bowtie2.sam")
  cmdline = '%s -x %s -1 %s -2 %s -U %s,%s -S %s -q --end-to-end --very-sensitive -p 12 --no-unal --phred%i --fr -I 0 -X 2000 --dovetail > %s_bowtie2.log 2> %s_bowtie2.err' % (BTBIN, db, fwd,rev,fwd_unpaired,rev_unpaired, samout, phredbase,prefix,prefix)

  print cmdline
  subprocess.call(cmdline, shell=True)
  

def postprocess_bowtie(prefix):
  # generate a sorted bam file from each samfile
  samname=os.path.join(ALDIR,prefix+"_bowtie2.sam")
  bamname_un = os.path.join(ALDIR,prefix+"_unsorted.bam")
  bamname=os.path.join(ALDIR,prefix+"_bowtie2_sorted")
  cmdline1 = "samtools view -bS %s > %s" % (samname,bamname_un)
  cmdline2 = "samtools sort %s %s" % (bamname_un,bamname)
  print cmdline1
  print cmdline2
  retcode1 = subprocess.call(cmdline1,shell=True)
  retcode2 = subprocess.call(cmdline2,shell=True)
  if retcode1 == 0 and retcode2 == 0:
    print "Safe to remove samfile %s" % prefix
    os.remove(samname)
    os.remove(bamname_un)
  else:
    print "*** Encountered an error while postprocessing %s" % prefix

# now we can actually run the samples

#pool1=Pool(processes=NPROC)
for samp in all_samples:
    #pool1.apply_async(preprocess_gz_file, samp)
    print "running on %s" % samp
    apply(preprocess_gz_file, [samp])

#pool1.close()
#pool1.join()

for samp in all_samples:
  run_bowtie(samp.outprefix,samp.phredbase,db=SEQ_DB)
  postprocess_bowtie(samp.outprefix)
