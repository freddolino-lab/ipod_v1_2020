#!/usr/bin/python

import matplotlib
matplotlib.use("Agg")
import numpy
import scipy.integrate
import scipy.stats
import os.path
import pylab

# calculate chip enrichment on each input file in my sequencing manifest
# we plot , and calculate the AUC for, the curve of normalized enrichment (against the maximum observed) vs cumulative distribution
# we assume that bootstrapping has already run

# useful constants
READDIR="bootstrap" ;# directory containing the trimmed reads and fastqc output

q_suffix = "_original.npy" ;# suffix for the file to be considered

def get_enr_stats(infile, outpng, outtxt):
    # read the given input file, and calculate the chip
    #  enrichment distribution
    # write an image to outpng and information to outtxt


    vals=numpy.load(infile)
    vals = vals/numpy.max(vals)

    quantiles = numpy.arange(100)

    scores = numpy.array( [scipy.stats.scoreatpercentile( vals,q) for q in quantiles] )

    pylab.figure()
    pylab.plot(quantiles,scores)
    pylab.savefig(outpng)

    # also calculate and write the auc
    auc=scipy.integrate.trapz(x=quantiles,y=scores)
    fout = open(outtxt,'w')
    fout.write("%f\n"% auc)
    fout.close()

if __name__ == '__main__':

    intab = open('read_manifest.txt')
    samp_prefixes = []
    for line in intab:
      if line[0] == '#':
          continue

      this_prefix = line.rstrip().split()[4]
      infile=os.path.join(READDIR, this_prefix + q_suffix)
      out_img = os.path.join(READDIR, this_prefix + '_chipenrich.png')
      out_txt = os.path.join(READDIR, this_prefix + '_chipenrich.txt')
      get_enr_stats(infile, out_img, out_txt)

    intab.close()
