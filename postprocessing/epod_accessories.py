# Accessory functions needed for EPOD calling

import numpy
import scipy
import scipy.stats
import scipy.signal
import bisect

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

def circular_range_bps(offsetvec_orig, datvec_orig, startbp, endbp, genomelength=None, guess_start=None, guess_end = None, return_inds = False):

  """
  Find all locations in a set of genomic coordinates and data that are within the given range, with appropriate wrapping assuming a circular chromosome

  offsetvec_orig should be the locations in bp, and datvec_orig the corresponding data of interest

  startbp and endbp are the 0-indexed beginning and end of the region to be considered

  This means that the number of elements returned is uncertain, since it
    depends on the spacing between entries

  We return the set of all entries completely contained with startbp
    and endbp, and the corresponding slice of offsetvec

  If genomelength is given, it is used as the total length of the genome
  Otherwise, we estimate it based on the last two entries

  guess_start and guess_end are guesses at where we should start looking,
    to avoid having to search through the entire array

  If return_inds is true, we return the indices passed to circular_range
    as well as the usual values
  """

  import bisect


  offsetvec = numpy.array(offsetvec_orig, copy=True)
  datvec = datvec_orig

  startbp=int(startbp)
  endbp=int(endbp)

  arrlen = len(datvec)
  if genomelength:
    lastbp = genomelength ;# last coordinate in the genome
  else:
    lastbp = offsetvec[-1] + (offsetvec[-1] - offsetvec[-2])

  if (startbp > endbp):
    tmpbp = endbp
    endbp=startbp
    startbp=tmpbp

    guesstmp = guess_start
    guess_start = guess_end
    guess_end = guesstmp

  if (startbp < 1):
    startbp += lastbp

  if endbp > offsetvec[-1]:
    endbp -= lastbp

  # limit the bounds of the search if we are able to
  searchleft_start = 0
  searchleft_end = len(offsetvec)
  searchright_start = 0
  searchright_end = len(offsetvec)
  width=2*abs(startbp-endbp)


  if (guess_start is not None) and (guess_start%arrlen > width) and (guess_start%arrlen < (arrlen-width)):
    searchleft_start = guess_start - width
    searchleft_end = guess_start + width
    #print "Using guessed start"
  if (guess_end is not None) and (guess_end%arrlen > width) and (guess_end%arrlen < (arrlen-width)):
    searchright_start = guess_end - width
    searchright_end = guess_end + width

  #print "Search positions are %i--%i %i--%i" % (searchleft_start, searchleft_end, searchright_start, searchright_end)

  startindex = bisect.bisect_left(offsetvec, startbp, searchleft_start, searchleft_end)
  endindex = bisect.bisect_right(offsetvec, endbp, searchright_start, searchright_end)

  if (startindex > endindex):
    offset_slice = numpy.concatenate((offsetvec[startindex:],offsetvec[:endindex]))
    data_slice = numpy.concatenate((datvec[startindex:],datvec[:endindex]))
  else:
    offset_slice = offsetvec[startindex:endindex]
    data_slice = datvec[startindex:endindex]

  #offset_slice = circular_range(offsetvec, startindex, endindex)
  #data_slice = circular_range(datvec, startindex, endindex)

  #print (startindex, endindex)

  if (return_inds):
    return (offset_slice, data_slice, startindex, endindex)
  else:
    return (offset_slice, data_slice)

def do_runningavg_opt(txtfile, outfile, width=49, genomelength=None):
  """
  Do a running average (mean) over a given width (in bp)

  This function **is** safe to jumps in bp numbering
  """

  if (width % 2 == 0):
    print "WARNING! Convolutions with an even width are not a good idea"

  halfwidth = width / 2

  instr1 = open(txtfile, 'r')
  instr2 = open(txtfile, 'r')
  ostr = open(outfile, 'w')

  # first figure out the total number of lines
  start = instr1.tell()
  nr = 0
  currline = instr1.readline()
  firstind, firstval = currline.split()
  firstind = int(firstind)
  finalind, finalval = (0,0)
  while(currline != ""):
    nr += 1
    if (currline != " "):
      lastind, lastval = (finalind, finalval)
      finalind, finalval = currline.split()
    currline = instr1.readline()


  if not genomelength:
    genomelength = int (finalind) + int(finalind) - int(lastind)
    print "Guessing genome length of %i" % genomelength

  instr1.seek(start)

  vallist = []; # list of the current values for the average
  offsetlist = []; # corresponding offsets

  # read the values for the wrapped end of the chromosome
  # Note that we search out the appropriate entry based on the difference in bp
  for i in range(nr - ( (width+1)/2)):
    instr1.readline()

  firstbp = firstind - halfwidth
  lastbp = firstind + halfwidth

  currline = instr1.readline()
  while (currline != ""):
    index, val = currline.split()
    index = int(index)
    if index > firstbp:
      vallist.append(float(val))
      offsetlist.append(int(index))
    currline = instr1.readline()

  instr1.seek(start)

  # use instr1 to load up the sets of values up to the end of the first
  #   moving average
  currline = instr1.readline()
  index, val = currline.split()
  index = int(index)
  while (index < lastbp):
    vallist.append(float(val))
    offsetlist.append(index)
    currline = instr1.readline()
    index, val = currline.split()
    index = int(index)

  #print "Calculating running average with %i elements" % len(vallist)
  #ostr.write("Current list of elements for vallist: \n")
  #for offset, elem in zip(offsetlist, vallist):
  #  ostr.write("%s %s\n" % (offset,elem))
  #ostr.write("----------\n")
  curravg = scipy.mean(vallist)
  currlen = len(vallist)

  # now start walking through the array and writing the correct values to ostr
  curr_center_line = instr2.readline()

  while (curr_center_line != ""):
    center_offset, center_val = curr_center_line.split()
    center_offset = int(center_offset)
    index, val = currline.split()
    index = int(index)
    val = float(val)

    goodrange = [center_offset - halfwidth, center_offset + halfwidth]
    #print "Centered on %s; range is %s - %s" % (curr_center_line, goodrange[0], goodrange[1])

    if goodrange[0] < 0:
      prevrange = [genomelength + goodrange[0], genomelength]
      goodrange[0] = 0
    elif goodrange[1] > genomelength:
      prevrange = [0, goodrange[1] - genomelength]
      goodrange[1] = genomelength
    else:
      prevrange = [-1,-1]

    #print "Acceptable ranges: %s %s" % (goodrange, prevrange)

    # find how big a slice of the list we need to cut
    removeind = -1
    for i in range(len(offsetlist)):
      curroffset = offsetlist[i]
      if (curroffset <= goodrange[1] and curroffset >= goodrange[0]):
        break
      if (curroffset <= prevrange[1] and curroffset >= prevrange[0]):
        break
      removeind = i
      
    #print removeind
    #print currlen
    #print vallist
    #print offsetlist
    #print "---"
    if (removeind >= 0):
      removeind += 1
      removedvals = vallist[:removeind]
      vallist = vallist[removeind:]
      offsetlist = offsetlist[removeind:]

      newlen = currlen - removeind
      if (newlen == 0):
        currlen = 0
        curravg = 0

      else:
        curravg = (currlen/float(newlen)) * curravg - numpy.sum(removedvals) / float(newlen)
        currlen=newlen

    # Add in new values as needed
    #if (index > genomelength):
    #  index -= genomelength

    if (index > genomelength):
        instr1.seek(start)
      
    #print "New index start: %i" % index
    while ( (index >= goodrange[0] and index <= goodrange[1]) or (index >= prevrange[0] and index <= prevrange[1])):
      vallist.append(val)
      newlen = currlen + 1
      curravg = curravg * (currlen / float(newlen)) + val / float(newlen)
      currlen = newlen
      offsetlist.append(index)
      currline = instr1.readline()
      #print "Current line: %s" % currline
      if (currline == ""):
        instr1.seek(start)
        currline = instr1.readline()
        #print "Wrapping to beginning"
      index, val = currline.split()
      index = int(index)
      val = float(val)

    #print "For center offset %i, lists are %s || %s" % (center_offset, vallist, offsetlist)
    #curravg = scipy.mean(vallist)
    ostr.write("%i %f\n" % ( int(center_offset), curravg ) )

    curr_center_line = instr2.readline()

  instr1.close()
  instr2.close()
  ostr.close()


def identify_epods_v3(epod_data, percentile_data, min_epod_length, epod_outfile, lpod=False, delta=25):
  """
  Look for epods over (100 - delta-th) percentile from percentile_data in epod_data

  Each data input should be a file (NOT a vector)
  We then look for all contiguous regions of length at least min_epod_length  
  A file (epod_outfile) containing the locations of epods will also be written

  If lpod is true, then we instead look for similar regions under 25th percentile occupancy

  All units are IN BASE PAIRS

  This version of the function tries to expand each epod symmetrically starting from local maxima , which is less likely
   to yield asymmetric epods than is the original identify_epods function

  In addition, drops in the occupancy trace below 0 will break the epod no matter what

  """

  offsets, epod_vec = numpy.loadtxt(epod_data, unpack=True)
  percentile_vec = numpy.loadtxt(percentile_data, usecols=(1,), unpack=True)
  epod_cutoff = scipy.stats.scoreatpercentile(percentile_vec, 100-delta)
  lpod_cutoff = scipy.stats.scoreatpercentile(percentile_vec, delta)
  percentile_vec = []

  # establish a guess for how many bp between each entry
  stride = offsets[1] - offsets[0]

  def epod_cmp(value):
    if (lpod):
      return (value < lpod_cutoff)
    else:
      return (value > epod_cutoff)

  if (lpod):
    print "Searching for regions at least %i bp long with a median z score below %f" % (min_epod_length, lpod_cutoff)
  else:
    print "Searching for regions at least %i bp long with a median z score above %f" % (min_epod_length, epod_cutoff)

  epod_pot_arr = numpy.zeros_like(epod_vec) - 1000 ;# contains the values of all 1024 bp windows that are potential epods
  # we follow a two-pass approach to find epods
  # first we go through the raw data and find all 1024 bp windows where the
  #   median value is at least percentile_vec
  # Each time we find such a window, we add to epod_pot_arr the score at that location
  #  this way we know the relative heights of various windows
  #  note that we add the plain score, and not the window median, to minimize ties

  offset = int(min_epod_length / 2)
  for i in range(len(offsets)):
    start = offsets[i]-offset
    end = offsets[i]+offset
    curr_median = scipy.median(circular_range_bps(offsets,epod_vec, start, end)[1])
    if (epod_cmp(curr_median)):
      epod_pot_arr[i] = epod_vec[i]

  # now, we find the BEST locations to start potential epods, and try expanding them in either direction

  #numpy.save('test_centers.npy',epod_pot_arr)
  epod_pot_centers = scipy.signal.argrelextrema(epod_pot_arr, numpy.greater_equal,mode='wrap')[0]
  epod_abovezero = numpy.argwhere(epod_pot_arr > -100)

  epod_centers = numpy.intersect1d(epod_pot_centers, epod_abovezero)

  epod_locs = []

  #print "DEBUG: Potential epod centers found at: %s" % epod_centers


  for center_i in epod_centers:
      # we start at just the centers, which are peaks in the occupancy trace, and make sure that we can expand to
      #  be large enough for an epod without hitting any zeroes
      epod_start = offsets[center_i] - 1
      epod_end = offsets[center_i] + 1

      expand_left = True
      expand_right = True

      while (expand_left or expand_right):
          # we try expanding this epod as far as we can

          if expand_left:
              # try expanding to the left
              # we break if either the window median drops too low, or the value at the position of interest drops below 0
              trial_median = scipy.median(circular_range_bps(offsets,epod_vec, epod_start - 1, epod_end)[1])
              new_value = circular_range_bps(offsets,epod_vec, epod_start - 1, epod_end)[1][0]
              if new_value < 0:
                  expand_left = False
              elif trial_median > epod_cutoff:
                  epod_start -= 1
              else:
                  expand_left = False

          if expand_right:
              # try expanding to the right
              trial_median = scipy.median(circular_range_bps(offsets,epod_vec, epod_start, epod_end+1)[1])
              new_value = circular_range_bps(offsets,epod_vec, epod_start - 1, epod_end)[1][-1]
              if new_value < 0:
                  expand_right = False
              elif trial_median > epod_cutoff:
                  epod_end += 1
              else:
                  expand_right = False
      
      #print "DEBUG: After expansion, ipod between %i and %i has median %f" % (epod_start, epod_end, scipy.median(circular_range_bps(offsets,epod_vec, epod_start, epod_end)[1]) )
      epod_locs.append( (int(epod_start), int(epod_end)) )
  
  # Take one shot at merging adjacent ipods
  i = 0
  while (i < (len(epod_locs) - 2)):
    j = i+1
    start1, end1 = epod_locs[i]
    start2, end2 = epod_locs[j]
    if (end1 > start2):
      #print "DEBUG: Trying to merge epods %s and %s" % (epod_locs[i], epod_locs[j])
      #print scipy.median(circular_range_bps(offsets,epod_vec, start1, end2)[1])
      if epod_cmp(scipy.median(circular_range_bps(offsets,epod_vec, start1, end2)[1])):
        #print "DEBUG: Merging epods %s and %s" % (epod_locs[i], epod_locs[j])
        epod_locs.pop(j)
        epod_locs[i] = (min(start1,start2),max(end1,end2))
        continue

    i += 1

  #print "DEBUG:" 
  #print epod_locs

  # write a file containing 1s at the positions involved in epods
  epod_loc_vec = scipy.zeros(len(epod_vec))
  for i in range(len(epod_vec)):
    iloc = offsets[i]
    for (start,end) in epod_locs:
      if numpy.abs(start-end) < min_epod_length:
          continue

      if (start > end):
        if (end < 0):
          end += len(epod_vec)
        else:
          start -= len(epod_vec)
      if (iloc>=start and iloc<= end):
        epod_loc_vec[i] = 1
        continue

  write_grfile(offsets, epod_loc_vec, epod_outfile)


  return epod_locs
