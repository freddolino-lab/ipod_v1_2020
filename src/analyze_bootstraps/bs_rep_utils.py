# general utility functions for working with bootstrap replicate matrices

# this has commonly used functions like normalizing the columns, or calculating CIs

import matplotlib
matplotlib.use("Agg")

import numpy
import scipy.stats
import scipy.interpolate
import ipod_utils
import sys
import pylab
import rpy2.robjects as robjects


def load_tall( infile ):
    # load a numpy array from a file, ensuring that it has at least two dimensions defined
    # if it does not, we add a second dimension (of length 1)
    # this is designed to coerce vectors into N x 1 matrices

    inmat = numpy.load(infile)

    if len(inmat.shape) == 1:
        inmat = inmat.reshape( (-1,1) )

    return inmat

    

def mean_norm_bs_mat(inmat,targetval=100.0, offset=0.25):
    # normalize a bootstrap matrix and return the normalized matrix
    # we assume that rows are positions and replicates are columns
    # each replicate is normalized as it would be by the normalize_mean_value_withoffset function
    #  from ipod_utils

    new_mat = numpy.array( inmat, copy=True)
    curr_means = numpy.mean(new_mat, axis=0)
    new_mat *= ( (targetval-offset) / curr_means)
    new_mat += offset
    return new_mat

def calc_bsmat_lograt(inmat1, inmat2):
    # calculate the log ratio between two matrices
    return numpy.log2(inmat1 / inmat2)

def get_ci_vecs(inmat):
    # give the element-wise top and bottoms of 95% cis for a bootstrap matrix
    # we return a tuple of numpy arrays with the bottom and top of that range

    ci_lo = scipy.stats.scoreatpercentile(inmat, 2.5, axis=1)
    ci_hi = scipy.stats.scoreatpercentile(inmat, 97.5, axis=1)

    return (ci_lo, ci_hi)

def get_elementwise_extrema( inmats ):
    # given a list of vectors, return vectors containing the elementwise minimum
    #  and maximum among the input values
    # 
    # We assume that each input matrix is an X by 1 matrix

    if len(inmats) == 0:
        raise(ValueError("Must give at least one vector to get_elementwise_extrema!"))

    if len(inmats) == 1:
        return (inmats[0],inmats[0])

    vlen,vheight = inmats[0].shape
    if vheight != 1:
        raise(ValueError("get_elementwise_extrema assumes N x 1 arrays"))

    for v in inmats[1:]:
        this_l, this_h = v.shape
        if (vlen != this_l):
            raise(ValueError("get_elementwise_extrema requires all arrays to have equal lengths"))

        if (this_h != 1):
            raise(ValueError("get_elementwise_extrema assumes N x 1 arrays"))

    # if we made it to this point, we can actually make one big matrix and then
    #  calculate the extrema

    full_mat = numpy.hstack( inmats )
    min_vals = numpy.min(full_mat, axis=1)
    max_vals = numpy.max(full_mat, axis=1)

    return( min_vals.reshape( (-1,1) ), max_vals.reshape( (-1,1) ) )

def calc_normalized_lograt_mats( inmat_1, inmat_2 ):
    # normalize the two input matrices, and then calculate and return their log ratio
    # we do NOT act on these matrices in place

    realmat_1 = numpy.array(inmat_1, copy=True)
    realmat_2 = numpy.array(inmat_2, copy=True)
    mean_norm_bs_mat(realmat_1)
    mean_norm_bs_mat(realmat_2)

    return calc_bsmat_lograt(realmat_1, realmat_2)

def do_percentile_normalization_bycol(inparr, percentile, targetval, offset):
    # apply normalization equivalent to that in
    # ipod_utils.normalize_percentile_value
    # each collumn is normalized separately
    # this acts on the input IN PLACE 

    for col in range(inparr.shape[1]):
        curr_quantile = scipy.stats.scoreatpercentile(inparr[:,col], percentile)
        inparr[:,col] *= ( (targetval - offset) / curr_quantile)
        inparr[:,col] += offset

def calc_lin_subtract_counts_log_vec( ipod_vec, chip_vec, soffset):
    """
    Helper function for calc_lin_subtract_counts_log_bycol

    Applies the chip-subtraction described there to a single row
    """

    sortord = numpy.argsort(chip_vec)
    ipod_forfit = ipod_vec[sortord]
    chip_forfit = chip_vec[sortord]

    goodflags = (chip_forfit > soffset)
    ipod_forfit = ipod_forfit[goodflags]
    chip_forfit = chip_forfit[goodflags]

    slope,intercept,rval,pval,stderr = scipy.stats.linregress( chip_forfit, ipod_forfit)
    zeroval = -1 * intercept / slope

    predvals = numpy.maximum(slope * chip_vec + intercept, numpy.zeros_like(chip_vec))

    newvals = ipod_vec - (predvals + soffset)
    return newvals





def calc_lin_subtract_counts_log_bycol( ipod_mat, chip_mat, soffset=0.0):
    """
    Equivalent to lin_subtract_counts_log from ipod_utils, applied column-wise

    do a linear fit to find the best mapping of chip_mat to ipod_mat, and then subtract from ipod_mat (note that we only fit things with an rna polymerase value over 0)
    soffset is added to all of the fitted values prior to subtraction
    """

    nrow,ncol = ipod_mat.shape
    if ipod_mat.shape != chip_mat.shape:
        raise(ValueError("IPOD and CHIP matrices must match in shape!"))

    output_mat = numpy.zeros_like(ipod_mat)

    for col_i in range(ncol):
        output_mat[:,col_i] = calc_lin_subtract_counts_log_vec( ipod_mat[:,col_i], chip_mat[:,col_i], soffset)

    return output_mat

def calc_v2x_chipsub_mats( ipod_mat, inp_mat, chip_mat, input_spline_trace):
    # calculate the chipsub_v2x_inplognormed sttatistic for a data set
    # each of the input matrices should be m x n, for m locations and n bootstrap replicates
    # we return a vector of n values corresponding to the chipsub value at each replicate, focusing on the
    # v2x_90per_lognorm method
    # we apply appropriate columnwise normalizations to a copy o the iarrays

    # constants for the v2x method
    V2X_PERC = 90.0
    V2X_TARGVAL = 500.0
    V2X_OFFSET = 0.25
    V2X_OFFSET_LATE = 0.0

    # first generate appropriately normalized inputs
    ipod_mat_90percnorm = numpy.array(ipod_mat, copy=True) / input_spline_trace
    # this is only for debugging
    ##numpy.save("ipod_test_splnormed.npy", ipod_mat_90percnorm)
    ##numpy.save("ipod_test_before_splnormed.npy", ipod_mat)
    # end debug

    inp_mat_90percnorm = numpy.array(inp_mat, copy=True) / input_spline_trace
    chip_mat_90percnorm = numpy.array(chip_mat, copy=True) / input_spline_trace

    do_percentile_normalization_bycol(ipod_mat_90percnorm, V2X_PERC, V2X_TARGVAL, V2X_OFFSET)
    do_percentile_normalization_bycol(inp_mat_90percnorm, V2X_PERC, V2X_TARGVAL, V2X_OFFSET)
    do_percentile_normalization_bycol(chip_mat_90percnorm, V2X_PERC, V2X_TARGVAL, V2X_OFFSET)

    # this is only for debugging
    ##numpy.save("ipod_test_med_normed.npy", ipod_mat_90percnorm)
    ##numpy.save("inp_test_med_normed.npy", inp_mat_90percnorm)
    ##numpy.save("chip_test_med_normed.npy", chip_mat_90percnorm)
    # end debug

    # do count subtraction and calculate log ratios vs input
    ipod_minus_chip = calc_lin_subtract_counts_log_bycol( ipod_mat_90percnorm, chip_mat_90percnorm, soffset=0.0)
    ## debug numpy.save("test_subtracted.npy", ipod_minus_chip)
    ipod_minus_chip[ipod_minus_chip < 0.25] = 0.25
    ipod_minus_chip_lognormed = calc_bsmat_lograt( ipod_minus_chip, inp_mat_90percnorm)
    ## debug numpy.save("test_sub_lognormed.npy", ipod_minus_chip_lognormed)

    # normalize and calculate log ratio for the chip subtracted file
    numpy.save("test_clamped.npy", ipod_minus_chip)
    do_percentile_normalization_bycol(ipod_minus_chip, V2X_PERC, V2X_TARGVAL, V2X_OFFSET_LATE)
    ## debug numpy.save("test_sub_percnormed.npy", ipod_minus_chip)
    ipod_minus_chip_90perclognorm = calc_bsmat_lograt( ipod_minus_chip, inp_mat_90percnorm)
    ## debug numpy.save("test_inp_90percnormed.npy", inp_mat_90percnorm)
    ## debug numpy.save("test_sub_percnormed_90perclognormed.npy", ipod_minus_chip_90perclognorm)

    return ipod_minus_chip_90perclognorm

def get_lin_sub_mat(mat_a, mat_b, plotfile=None, minperc=98.0):
    # apply a linear regression to subtract predictions based on mat_b from mat_a
    # we apply a constrained linear regression which must pass through the origin, acting only on points with mat_b> minval, and
    #  then use a slope that keeps 95% of those high-ChIP samples below the line
    # minperc is the minimum percentile of the chip data (mat_b) to use

    minval = scipy.stats.scoreatpercentile(mat_b, minperc)
    print "Minval is %f" % minval
    out_mat = numpy.zeros_like(mat_a)

    if plotfile is not None:
        pylab.figure()
        pylab.rcParams.update({'font.size': 22})
        pylab.hexbin(mat_b, mat_a, bins='log', cmap=pylab.get_cmap("Blues"))
        pylab.xlabel('Log$_2$ ChIP vs. input')
        pylab.ylabel('Log$_2$ IPOD vs. input')
        pylab.savefig(plotfile+"_noline.png")
        #pylab.plot(col_b[numpy.argsort(col_b)], pred_for_plot, 'r-', linewidth=2.5)

    for col_i in range(mat_a.shape[1]):

        col_a = mat_a[:,col_i]
        col_b = mat_b[:,col_i]

        fitslope =do_lin_fit(col_b, col_a, minval)
        print "Using slope of %f" % fitslope
        predvals = do_lin_interpolation_array(col_b,fitslope, 0)

        predvals = numpy.fmax(predvals, 0.0)

        newvals = col_a - predvals

        print col_b.shape
        print predvals.shape
        print "***"
        if (plotfile is not None) and (col_i == 0):
            pylab.plot(col_b[numpy.argsort(col_b)], predvals[numpy.argsort(col_b)], 'r--', linewidth=3)
            pylab.savefig(plotfile + "_main.png")
            pylab.figure()
            pylab.hexbin( col_a, newvals, bins='log', cmap=pylab.get_cmap("Blues"))
            pylab.colorbar()
            pylab.xlabel('Log$_2$ IPOD vs. input')
            pylab.ylabel('ChIP-subtracted IPOD-HR')
            pylab.savefig(plotfile + "_sub.png")
            pylab.figure()
            pylab.hexbin( col_b, newvals, bins='log', cmap=pylab.get_cmap("Blues"))
            pylab.colorbar()
            pylab.xlabel('Log$_2$ ChIP vs. input')
            pylab.ylabel('ChIP-subtracted IPOD-HR')
            pylab.savefig(plotfile + "_chipsub_lin.png")

        #if len(newvals.shape) == 1:
        #    newvals = newvals.reshape(-1,1)

        out_mat[:,col_i] = newvals


    return out_mat

def get_rzscores_bycol(inp_mat):
    # return the robust z score normalized version of a matrix, acting separately on each column

    ncol = inp_mat.shape[1]
    outmat = numpy.zeros_like(inp_mat)

    for i in range(ncol):
        this_col = inp_mat[:,i]
        this_mad = 1.4826 * numpy.median( numpy.abs( this_col - numpy.median(this_col) ) )
        zscores = (this_col - numpy.median(this_col)) / this_mad
        outmat[:,i] = zscores

    return outmat


def calc_v6_chipsub_mats( ipod_vs_inp_mat, chip_vs_inp_mat, plotfile=None):
    # do v6 chip subtraction following the method in
    #  my chip subtraction script
    # each of the input matrices should be m x n, for m locations and n bootstrap replicates
    # we return a matrix of mxn values corresponding to the chipsub value at each replicate, containing the rscores
    # recall that we want as input the ipod_vs_inp and chip_vs_inp log ratios

    # first calculte the count subtraction
    # this uses a linear fit that keeps the majority of the high-ChIP data below the line

    sub_mat = get_lin_sub_mat( ipod_vs_inp_mat, chip_vs_inp_mat, plotfile )


    # following line is for DEBUG ONLY
    numpy.save('sub_mat.npy',sub_mat)

    # now convert to robust z-scores and return
    zscore_mat = get_rzscores_bycol( sub_mat )

    return zscore_mat

def collapse_rep_matrices_mean_norm( inp_files ):
    # collapse two or more input matrices into a single output matrix, averaging their entries
    # we mean-normalize the matrices prior to combining them
    # inp_locs should contain the locations of all of the probes
    # then return the normalized, merged matrix

    orig_vals = [numpy.load(f) for f in inp_files]
    for i in range(len(orig_vals)):
        if len(orig_vals[i].shape) == 1:
         orig_vals[i] = orig_vals[i].reshape( (-1,1) )

    orig_normed = [ mean_norm_bs_mat(i,offset=0) for i in orig_vals ]
    normed_vals = reduce( lambda x,y: x+y, orig_normed[1:], orig_normed[0]) / float(len(orig_normed))
    return normed_vals

def collapse_rep_matrices_get_extrema( inp_files ):
    # apply mean normalization to a set of matrices, and then return the elementwise extrme
    # we assume that each input matrix is N x 1, as in get_elementwise_extrema
    # at present we do NOT re-normalize the final matrix, although I need to consider this carefully

    orig_vals = [numpy.load(f).reshape( (-1,1) ) for f in inp_files]
    orig_normed = [ mean_norm_bs_mat(i,offset=0) for i in orig_vals ]
    return get_elementwise_extrema(orig_normed)

def get_spl_norm( x_locs, inp_files, resolution, genome_length = 4641652, oriCloc=3923883 ):
    # generate a spline-smoothed version over a given set of input files
    #
    #  the genome_length and oriCloc arguments default to suitable values for U00096.3
    #  x_locs should be the locations, in bp, of all data points
    #  resolution is the spacing, in bp, between data points

    # first just collapse the replicates and normalize as we need to
    orig_vals = [numpy.load(f) for f in inp_files]
    orig_normed = [ mean_norm_bs_mat(i,offset=0) for i in orig_vals ]
    sum_vals = reduce( lambda x,y: x+y, orig_normed[1:], orig_normed[0]) / float(len(orig_normed))

    # now actually do the spline fit and return those values
    locs_recentered = (x_locs - oriCloc) % genome_length
    sortorder = numpy.argsort(locs_recentered)
    locs_mod = locs_recentered[sortorder]
    vals_mod = sum_vals[sortorder]
    locs_orig_mod = x_locs[sortorder]
    resort_order = numpy.argsort(locs_orig_mod)

    locs_aug = numpy.append(locs_mod, genome_length + locs_mod[0])
    vals_aug = numpy.append(vals_mod, vals_mod[0])
                                                                                            
    knots=[genome_length / 4, genome_length/2, 3 * genome_length / 4]  
    allspl = scipy.interpolate.splrep(locs_aug, vals_aug, k = 3, t=knots, task=-1, w=scipy.ones(len(locs_aug)), per=1, full_output=1) 
    myspl=allspl[0]
    if allspl[2] > 0:
        raise("ERROR FITTING SPLINE") 

    splvals = (scipy.interpolate.splev(locs_mod, myspl))[resort_order]
    splvals /= numpy.mean(splvals)
    splvals = splvals.reshape( (-1,1) )

    # the plotting here is for debugging purposes only
    #pylab.figure()
    #pylab.plot(x_locs, sum_vals, 'b-')
    #pylab.plot(x_locs, splvals, 'r-')
    #pylab.savefig('spl.png')
    
    return splvals

def get_all_comb_vals( function, this_args, other_arrs ):
    """
    Recursive function to enumerate all possible combinations of function arguments

    If other_arrs has length 1, return a list containing the
    result of function() by combining the arguments in this_args with each argument in the list given by other_arrs
    Otherwise, prepare a list consisting of all possible combinations of values in other_arrs[0] with downstream arguments

    This is intended mainly as a helper function for the replicate-level combination code below
    """

    if len(other_arrs) < 1:
        raise(ValueError("Don't know how to handle an empty array in get_all_comb_vals"))

    if len(other_arrs) == 1:
        retvals = []
        for val in other_arrs[0]:
            all_args = this_args + tuple([val])
            retvals.append(apply(function,all_args))

        return retvals

    retvals = []
    for val in other_arrs[0]:
        new_args = this_args + tuple([val])
        retvals += get_all_comb_vals(function, new_args, other_arrs[1:])

    return retvals

def get_replevel_extrema( function, dat_arrs):
    """
    return the smallest/largest possible values of function, in vector form

    evaluates function for each possible combination of values in 
    dat_arrs, and at each element gets the highest and lowest possible values

    function should take a number of aguments eual to the number
    of elements in dat_arrs, in the order specified in dat_arrs

    dat_arrs should be a list of lists; each fo the sublists should have one data vector per replicate 9presumably these are separate ipod, inp, and chip replicates, although we can handle a more general case here
    """

    # first we use a helper function to recursively generate all of the possible result vectors for different combinations of data arrays

    all_vals_list= get_all_comb_vals( function, (), dat_arrs)

    # now find the position-wise extrema
    all_vals_mat = numpy.array(all_vals_list)
    print "getting mat shape"
    print all_vals_mat.shape
    min_vals = numpy.amin(all_vals_mat, axis=0)
    max_vals = numpy.amax(all_vals_mat, axis=0)
    return (min_vals, max_vals)

def get_replevel_allpos( function, dat_arrs):
    """
    return an array of the results of function for all possible
      combinations of dat_arrs

    This is very similar to get_replevel_extrema, but returns all of the values
    """

    return get_all_comb_vals( function, (), dat_arrs)


def calc_normed_lograt(file_1, file_2, spl_vec):
    """
    calculate the normalized log ratio between two files and return it

    This functions very similarly to calc_ipod_inp_lograts, but
    we assume that each input file is a SINGLE file, and
    we subsequently return the resulting vector
    """

    vec_1 = numpy.load(file_1).reshape((-1,1)) / spl_vec
    vec_2 = numpy.load(file_2).reshape((-1,1)) / spl_vec

    vec_1_normed = mean_norm_bs_mat(vec_1)
    vec_2_normed = mean_norm_bs_mat(vec_2)

    this_lograt = calc_bsmat_lograt(vec_1_normed, vec_2_normed)

    return this_lograt

def calc_signed_log10p(vals):
    # given a set of values that can be considered z-score like,
    #   return the -log10 p value for each score in a standard 
    #   normal distribution

    pvals_raw = scipy.stats.norm.sf( vals)
    pvals_log = numpy.log10(pvals_raw)
    pvals_signed = (-1 * pvals_log)
    return pvals_signed

def do_lin_fit(xdat, ydat, minval):
  """
  do a linear regression on the selected data where xdat>minval, 
  and find and return the minimal slope for which 95% of the fitted data are below the line
  """

  print "Doing fit with minval %f" % (minval,)

  goodflags = (xdat > minval)

  robjects.r.assign('x', robjects.FloatVector(xdat[goodflags]))
  robjects.r.assign('y', robjects.FloatVector(ydat[goodflags]))
  robjects.r('my.linfit = lm(y ~ x+0)')
  robjects.r('slope.init=coef(my.linfit)')
  robjects.r('print(slope.init)')
  robjects.r('while ( sum( y > (slope.init * x) ) > (0.05 * length( y ) ) ) { slope.init = slope.init + 0.00001  } ')
  v=robjects.r('slope.init')

  print(v)
  return numpy.asarray(v)

def do_lin_interpolation(xval, slope, padding):
  """
  Return the predicted value for xval based on a linear regression with slope slope and intercept 0, plus padding
  We return 0 for negative values of x
  """

  return ((xval * slope + padding) > 0) * (xval * slope + padding)

def do_lin_interpolation_array(xvals, slope, padding):
  """
  Return the predicted value for each point in xvals based on a linear regression with slope slope and intercept 0, plus padding
  We return 0 for negative values of x
  """

  return ((xvals * slope + padding) > 0) * (xvals * slope + padding)

