#!/usr/bin/python

# apply the cwt-based peak caller from scipy.signal to a gr file

# usage:
#  python call_peaks_cwt.py ingrfile outgrfile threshold padsize
#  here threshold is the value to be passed to the peak caller's min_snr argument
#  padsize is the number of bins in either direction of each peak center to be added to the peak

import sys
import os.path
sys.path.append( os.path.join(os.path.realpath( __file__ ), "analyze_bootstraps" ) )
import numpy
import scipy.signal

def read_grfile(filename, skiprows=0):
  """
  Read data from a .gr file into a numpy array

  In the process we transpose the array, so that it has a[0] the indices and a[1] the values
  """

  locs = numpy.loadtxt(filename, dtype='int', usecols=(0,), skiprows=skiprows)
  vals = numpy.loadtxt(filename, usecols=(1,), skiprows=skiprows)
  return (locs,vals)


def write_grfile(indices, data, filename, header=None):
  """
  Write a gr file given the list of sequence positions and corresponding scalars
  """

  ostr = open(filename, 'w')
  if header is not None:
    ostr.write(header)
  for ind, dat in zip(indices, data):
    ostr.write("%i %f\n" % (ind, dat))
  ostr.close()


infile,outfile,snr_thresh, pad_size = sys.argv[1:]
snr_threshold=float(snr_thresh)
half_windowsize=int(pad_size)

locs,vals=read_grfile(infile)

peakwidths = numpy.arange(5,25)


peakinds = scipy.signal.find_peaks_cwt(vals, widths=peakwidths, min_snr = snr_threshold)
peakflags = numpy.zeros_like(vals)

for ind in peakinds:
    peakflags[ max((ind - half_windowsize),0) : min((ind + half_windowsize+1), len(peakflags)) ] = 1

write_grfile(locs,peakflags,outfile + ".gr")

