#!/usr/bin/python

# Tools for working with peak definitions

class PeakList:
  """
  Class for storing information on a set of peak locations

  Really all this means is that one has a collection of ranges in the genome, possibly with associated scalar data
  This is useful for looking at genes, ipod/chip peaks, etc.

  The primary focus of this implementation is for ipod peaks

  Note that we use closed intervals on both sides and 0 indexing for the internal representation
  """

  def __init__(self):
    self.peakstarts = []
    self.peakends = []
    self.peakvals = []
    self.numpeaks = 0

  def add_peak(self, start, end, value):

    self.peakstarts.append(start)
    self.peakends.append(end)
    self.peakvals.append(value)
    self.numpeaks += 1

  def get_peaks(self,withdat=False):

    if withdat:
      return zip(self.peakstarts, self.peakends, self.peakvals)
    else:
      return zip(self.peakstarts, self.peakends)

  def diff_peaks(self,other):
    """
    Return a new PeakList object containing  peaks in self that are not in other

    Order is maintained, and sorting is not assumed
    """

    newlist = PeakList()

    for (start,end,val) in self.get_peaks(True):
      is_unique = True
      for (ostart, oend) in other.get_peaks():
        if (check_peak_overlap(start,end,ostart,oend)):
          is_unique = False
          break

      if (is_unique):
        newlist.add_peak(start,end,val)

    return newlist

def read_flaggrfile(infile):
  """
  Read a gr file containing a list of 0 (no peak) or 1 (peak) flags and convert to a PeakList (returned)
  """

  instr = open(infile)
  
  is_peak = False
  currstart = 0
  currend = 0

  allpeaks = PeakList()

  prev_offset = 0

  for line in instr:
    linearr = line.rstrip().split()
    offset = int( round(float(linearr[0])) + 0.2)
    flag = int( round(float(linearr[1])) + 0.2)

    if is_peak:

      if flag > 0:
        prev_offset = offset
        continue

      else:
        currend = prev_offset
        allpeaks.add_peak(currstart, prev_offset, 0)
        prev_offset = offset
        is_peak = False
        continue

    else:

      if flag > 0:
        currstart = offset
        is_peak = True
        prev_offset = offset
        continue

      else:
        prev_offset = offset
        continue

    raise Error("We should never get here")

  # Handle any peaks wrapping around the origin
  if is_peak:
    print "Warning: Truncating peak at end of file"
    allpeaks.add_peak(currstart, prev_offset, 0)

  instr.close()
  print allpeaks.numpeaks
  return allpeaks

def check_peak_overlap(start1,end1,start2,end2):
  #Return True iff the peaks overlap
  if ( (start1 >= start2 and start1 <= end2) or (end1 >= start2 and end1 <= end2) or (start2 >= start1 and start2 <= end1) or (end2 >= start1 and end2 <= end1)):
    return True

  return False

