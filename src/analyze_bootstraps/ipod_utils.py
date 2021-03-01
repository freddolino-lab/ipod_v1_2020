#!/usr/bin/python

# A collection of functions useful in working with data from IPOD

import numpy

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

def read_grfile(filename, skiprows=0):
  """
  Read data from a .gr file into a numpy array

  In the process we transpose the array, so that it has a[0] the indices and a[1] the values
  """

  locs = numpy.loadtxt(filename, dtype='int', usecols=(0,), skiprows=skiprows)
  vals = numpy.loadtxt(filename, usecols=(1,), skiprows=skiprows)
  return (locs,vals)

