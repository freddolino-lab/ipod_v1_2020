#!/usr/bin/python

# do all of the preprocessing needed to apply Mike's bootstrapping script to my data sets
# we run bootstrapping to generate 200 samples at 5 bp resolution, and then write the results

# we require samtools v 0.1.19 for this to work. Newer versions may not function properly

import sys
import tempfile
import subprocess
import os.path

BINDIR=os.path.dirname( os.path.realpath( __file__ ) )

PARSE_CMD = "python %s/bootstrap_sam_file.py parse --paired %s %s"
SAMPLE_CMD = "python %s/bootstrap_sam_file.py  sample %s.npy %s.npy 4641652 --num_samples 200 --resolution 5"
ORIG_CMD = "python %s/bootstrap_sam_file.py  sample %s.npy %s.npy 4641652 --identity --resolution 5"

infile=sys.argv[1]
outprefix=sys.argv[2]
nthreads=int(sys.argv[3])


# first make a temporary, sorted sam file
tmp_file = tempfile.NamedTemporaryFile()
filter_cmd = "samtools sort -n -o -@ %i %s %s_tmp | samtools view -f 3 -F 2308 -q 30 - | sed s/\"#0\/4\t\"/\"#0\t\"/g> %s" % (nthreads, infile, infile, tmp_file.name)
#filter_cmd = "samtools sort -n -o -@ %i %s %s_tmp | samtools view -f 3 -F 2308 -q 30 - > %s" % (nthreads, infile, infile, 'test.bam')
print filter_cmd
subprocess.call(filter_cmd, shell=True)

# then actually run mike's command
print (PARSE_CMD % (BINDIR,tmp_file.name, outprefix+"_parsed"))
subprocess.call(PARSE_CMD % (BINDIR,tmp_file.name, outprefix+"_parsed"),shell=True)

# now do the sampling
print (SAMPLE_CMD % (BINDIR,outprefix+"_parsed", outprefix+"_resampled"))
subprocess.call(SAMPLE_CMD % (BINDIR,outprefix+"_parsed", outprefix+"_resampled"),shell=True)

# also generate a file with the original occupancies
print (ORIG_CMD % (BINDIR,outprefix+"_parsed", outprefix+"_original"))
subprocess.call(ORIG_CMD % (BINDIR,outprefix+"_parsed", outprefix+"_original"),shell=True)

# then clean up
tmp_file.close()
