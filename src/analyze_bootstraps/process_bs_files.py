#!/usr/bin/python

import numpy
from bs_rep_utils import *
import ipod_utils
import toml
import os.path
import sys


# this is obviously specific to E coli MG1655
G_LOCS = numpy.arange(0,4641652,5)

# functions to apply my usual ipod data processing to a bunch of bootstrap samples
# in the end this should allow me to collapse all of the replicates from ipod, chip, and
#  inp samples to my normal set of five summary files:
#  -ipod_vs_inp log2ratio signal
#  -chip_vs_inp log2ratio signal
#  -v2x chip subtracted log2ratio signal
#  -v6 rzscore signal
#  -v6 signed log p values
#  For each of these, I generate the observed value plus a 95% CI
# we assume that files called PREFIX_qnorm.npy exist
#  and have the quantile normalized occupancy data

def calc_ipod_inp_lograts( ipod_orig_files, ipod_bs_files, inp_orig_files, inp_bs_files, outprefix, spl_vec):
    # calculate observed and 95% ci values for the log ratio between ipod and inp samples
    #   the data are spline-corrected and mean-normalized prior to ratio calculation
    #  
    # the spline correction arises by dividing by spl_vec, which should be a spline-smoothed
    #   version of the occupancy from the input file
    # 
    # Note that this could be used equally well to calculate chip_vs_inp
    # We generate the following outputs, all prefixed with outprefix #   _lograt.gr -- the log ratio itself
    #   _lograt_95cilo.gr -- the bottom of a 95% ci on the log ratio
    #   _lograt_95cihi.gr -- the top of a 95% ci on the log ratio
    #   _lograt.npy -- log ratio of the observed values
    #   _lograt_bs.npy -- log ratios for each of the bootstrap samplings
    #   _lograt_replo.gr -- log ratio that is the lowest possible given replicate values
    #   _lograt_rephi.gr -- log ratio that is the highest possible given replicate values
    #   _lograt_replo.npy -- log ratio that is the lowest possible given replicate values
    #   _lograt_rephi.npy -- log ratio that is the highest possible given replicate values


    ipod_orig_vecs = [numpy.load(f).reshape((-1,1)) / spl_vec for f in ipod_orig_files]
    inp_orig_vecs = [numpy.load(f).reshape((-1,1)) / spl_vec for f in inp_orig_files]

    ipod_bs_mat_list = [numpy.load( f ) / spl_vec for f in ipod_bs_files ]
    inp_bs_mat_list = [numpy.load( f ) / spl_vec for f in inp_bs_files ]

    # normalize all of the matrices appropriately
    ipod_orig_normed = [ mean_norm_bs_mat(i) for i in ipod_orig_vecs ]
    inp_orig_normed = [ mean_norm_bs_mat(i) for i in inp_orig_vecs ]

    ipod_bs_normed = [ mean_norm_bs_mat(i) for i in ipod_bs_mat_list ]
    inp_bs_normed = [ mean_norm_bs_mat(i) for i in inp_bs_mat_list ]


    # collapse each data set across replicates

    ipod_orig_vals = reduce( lambda x,y: x+y, ipod_orig_normed[1:], ipod_orig_normed[0]) / float(len(ipod_orig_normed))
    inp_orig_vals = reduce( lambda x,y: x+y, inp_orig_normed[1:], inp_orig_normed[0]) / float(len(inp_orig_normed))

    ipod_bs_vals = reduce( lambda x,y: x+y, ipod_bs_normed[1:], ipod_bs_normed[0]) / float(len(ipod_bs_normed))
    inp_bs_vals = reduce( lambda x,y: x+y, inp_bs_normed[1:], inp_bs_normed[0]) / float(len(inp_bs_normed))

    ipod_rep_vals_min,ipod_rep_vals_max = get_elementwise_extrema( ipod_orig_normed )
    inp_rep_vals_min, inp_rep_vals_max = get_elementwise_extrema( inp_orig_normed )

    # calculate the requested log ratios, for the original values and alternates

    orig_lograt = calc_bsmat_lograt(ipod_orig_vals,inp_orig_vals )
    bs_lograt = calc_bsmat_lograt(ipod_bs_vals,inp_bs_vals )

    ci_lo, ci_hi = get_ci_vecs(bs_lograt)

    min_logmat = calc_bsmat_lograt( ipod_rep_vals_min, inp_rep_vals_max )
    max_logmat = calc_bsmat_lograt( ipod_rep_vals_max, inp_rep_vals_min )

    # now write the needed output
    ipod_utils.write_grfile(G_LOCS, orig_lograt, "%s_lograt.gr" % outprefix)
    ipod_utils.write_grfile(G_LOCS, ipod_orig_normed[0], "%s_DEBUG_ipod_orig_normed.gr" % outprefix)
    ipod_utils.write_grfile(G_LOCS, inp_orig_normed[0], "%s_DEBUG_inp_orig_normed.gr" % outprefix)
    ipod_utils.write_grfile(G_LOCS, ci_lo, "%s_lograt_95cilo.gr" % outprefix)
    ipod_utils.write_grfile(G_LOCS, ci_hi, "%s_lograt_95cihi.gr" % outprefix)
    ipod_utils.write_grfile(G_LOCS, min_logmat, "%s_lograt_replo.gr" % outprefix)
    ipod_utils.write_grfile(G_LOCS, max_logmat, "%s_lograt_rephi.gr" % outprefix)

    numpy.save("%s_lograt_bs.npy" % outprefix, bs_lograt)
    numpy.save("%s_lograt.npy" % outprefix, orig_lograt)
    numpy.save("%s_lograt_replo.npy" % outprefix, min_logmat)
    numpy.save("%s_lograt_rephi.npy" % outprefix, max_logmat)


# now actually implement the full processing of the bootstrap replicates

if __name__ == "__main__":


    conf_file = sys.argv[1]
    conf_dict = toml.load(conf_file)

    # figure out some global parameters
    bs_suffix = conf_dict['general']['bs_suffix']
    bs_dir = conf_dict['general']['bs_dirname']
    orig_suffix = conf_dict['general']['orig_suffix']
    out_prefix = os.path.join(conf_dict['general']['output_path'],conf_dict['general']['out_prefix'])

    # make missing path if needed
    if not(os.path.isdir(conf_dict['general']['output_path'])):
        os.mkdir(conf_dict['general']['output_path'])

    # get the spline normalization from the input data only
    #   we will use this repeatedly in the analysis below

    inp_dir = conf_dict['inp']['directory']
    inp_prefixes = conf_dict['inp']['sample_names']

    orig_pref,orig_suff = os.path.splitext(orig_suffix)
    orig_name = orig_pref + '_qnorm' + orig_suff

    inp_orig_files = [ os.path.join( inp_dir, bs_dir, x + orig_name) for x in inp_prefixes]
    input_spline_trace = get_spl_norm(G_LOCS, inp_orig_files, 5)


    # generate the collapsed ipod_vs_inp log ratios
    ipod_dir = conf_dict['ipod']['directory']
    ipod_prefixes = conf_dict['ipod']['sample_names']

    bs_pref,bs_suff = os.path.splitext(bs_suffix)
    bs_name = bs_pref + '_qnorm' + bs_suff

    ipod_orig_files = [ os.path.join( ipod_dir, bs_dir, x + orig_name) for x in ipod_prefixes]
    ipod_bs_files = [ os.path.join( ipod_dir, bs_dir, x + bs_name) for x in ipod_prefixes]

    inp_bs_files = [ os.path.join( inp_dir, bs_dir, x + bs_name) for x in inp_prefixes]

    calc_ipod_inp_lograts( ipod_orig_files, ipod_bs_files, inp_orig_files, inp_bs_files, out_prefix + "_ipod_vs_inp", input_spline_trace)

    # next calculate the collapsed chip_vs_inp log ratios
    chip_dir = conf_dict['chip']['directory']
    chip_prefixes = conf_dict['chip']['sample_names']

    chip_orig_files = [ os.path.join( chip_dir, bs_dir, x + orig_name) for x in chip_prefixes]
    chip_bs_files = [ os.path.join( chip_dir, bs_dir, x + bs_name) for x in chip_prefixes]

    calc_ipod_inp_lograts( chip_orig_files, chip_bs_files, inp_orig_files, inp_bs_files, out_prefix + "_chip_vs_inp", input_spline_trace)

    # do the v2x chip subtraction, building off of the results of the prior steps
    ## to get to this, we have to first collapse all of our replicates,
    ## and the files containing the original scores

    ipod_orig_merged  = collapse_rep_matrices_mean_norm( ipod_orig_files )
    inp_orig_merged  = collapse_rep_matrices_mean_norm( inp_orig_files )
    chip_orig_merged  = collapse_rep_matrices_mean_norm( chip_orig_files )

    ipod_bs_merged  = collapse_rep_matrices_mean_norm( ipod_bs_files )
    inp_bs_merged  = collapse_rep_matrices_mean_norm( inp_bs_files )
    chip_bs_merged  = collapse_rep_matrices_mean_norm( chip_bs_files )

    ipod_orig_min, ipod_orig_max = collapse_rep_matrices_get_extrema( ipod_orig_files )
    inp_orig_min, inp_orig_max = collapse_rep_matrices_get_extrema( inp_orig_files )
    chip_orig_min, chip_orig_max = collapse_rep_matrices_get_extrema( chip_orig_files )

    ## now actually do the subtraction on both the original and bootstrap matrix
    # DEBUG
    #numpy.save('ipod_orig_min.npy',ipod_orig_min)
    #numpy.save('ipod_orig_max.npy', ipod_orig_max)
    #numpy.save('inp_orig_min.npy',inp_orig_min)
    #numpy.save('inp_orig_max.npy',inp_orig_max)
    #numpy.save('chip_orig_min.npy',chip_orig_min)
    #numpy.save('chip_orig_max.npy', chip_orig_max)
    # /DEBug

    v2x_mat_bs = calc_v2x_chipsub_mats(ipod_bs_merged, inp_bs_merged, chip_bs_merged,input_spline_trace )
    #= calc_v2x_chipsub_mats( ipod_orig_min, inp_orig_max, chip_orig_max,input_spline_trace)
    #v2x_mat_rep_max = calc_v2x_chipsub_mats( ipod_orig_max, inp_orig_min, chip_orig_min,input_spline_trace)
    v2x_mat_orig = calc_v2x_chipsub_mats(ipod_orig_merged, inp_orig_merged, chip_orig_merged,input_spline_trace )

    # also run to get the most extreme values possible with these replicates
    ipod_normed_mats = [ mean_norm_bs_mat(numpy.load(f).reshape( (-1,1) ), offset=0) for f in ipod_orig_files ]
    inp_normed_mats = [ mean_norm_bs_mat(numpy.load(f).reshape( (-1,1) ), offset=0) for f in inp_orig_files ]
    chip_normed_mats = [ mean_norm_bs_mat(numpy.load(f).reshape( (-1,1) ), offset=0) for f in chip_orig_files ]
    v2x_mat_rep_min,v2x_mat_rep_max = get_replevel_extrema( calc_v2x_chipsub_mats, [ipod_normed_mats, inp_normed_mats, chip_normed_mats, [input_spline_trace] ] )

    ## write appropriate v2x output files
    numpy.save("%s_v2x_chipsub_orig.npy" % out_prefix, v2x_mat_orig)
    numpy.save("%s_v2x_chipsub_boot.npy" % out_prefix, v2x_mat_bs)
    numpy.save("%s_v2x_chipsub_rep_min.npy" % out_prefix, v2x_mat_rep_min)
    numpy.save("%s_v2x_chipsub_rep_max.npy" % out_prefix, v2x_mat_rep_max)

    ci_lo_v2x, ci_hi_v2x = get_ci_vecs(v2x_mat_bs)
    ipod_utils.write_grfile(G_LOCS, v2x_mat_orig, "%s_v2x_chipsub.gr" % out_prefix)
    ipod_utils.write_grfile(G_LOCS, ci_lo_v2x, "%s_v2x_chipsub_95cilo.gr" % out_prefix)
    ipod_utils.write_grfile(G_LOCS, ci_hi_v2x, "%s_v2x_chipsub_95cihi.gr" % out_prefix)
    ipod_utils.write_grfile(G_LOCS, v2x_mat_rep_min, "%s_v2x_chipsub_rep_min.gr" % out_prefix)
    ipod_utils.write_grfile(G_LOCS, v2x_mat_rep_max, "%s_v2x_chipsub_rep_max.gr" % out_prefix)

    # Now run through the full v6 analysis
    ## this is supposed to build off of the ipod and chip log ratios
    ## first we load the matrices from those input files

    ipod_vs_inp_orig = load_tall(out_prefix + "_ipod_vs_inp_lograt.npy")
    ipod_vs_inp_replo = load_tall(out_prefix + "_ipod_vs_inp_lograt_replo.npy")
    ipod_vs_inp_rephi = load_tall(out_prefix + "_ipod_vs_inp_lograt_rephi.npy")
    ipod_vs_inp_boot = load_tall(out_prefix + "_ipod_vs_inp_lograt_bs.npy")
    #debug numpy.save('debug_ipod_boot_mat.npy', ipod_vs_inp_boot)

    chip_vs_inp_orig = load_tall(out_prefix + "_chip_vs_inp_lograt.npy")
    chip_vs_inp_replo = load_tall(out_prefix + "_chip_vs_inp_lograt_replo.npy")
    chip_vs_inp_rephi = load_tall(out_prefix + "_chip_vs_inp_lograt_rephi.npy")
    chip_vs_inp_boot = load_tall(out_prefix + "_chip_vs_inp_lograt_bs.npy")
    #debug numpy.save('debug_chip_boot_mat.npy', chip_vs_inp_boot)

    v6_mat_orig = calc_v6_chipsub_mats(ipod_vs_inp_orig, chip_vs_inp_orig, plotfile = out_prefix + "_v6_linfit_sub.png") 
    v6_mat_boot = calc_v6_chipsub_mats(ipod_vs_inp_boot, chip_vs_inp_boot) 
    #debug numpy.save('debug_v6_boot_mat.npy', v6_mat_boot)

    # also have to do a little extra work to get the most extreme rep values

    all_ipod_vs_inp_vals = get_replevel_allpos( calc_normed_lograt, [ ipod_orig_files, inp_orig_files, [input_spline_trace]] )
    all_chip_vs_inp_vals = get_replevel_allpos( calc_normed_lograt, [ chip_orig_files, inp_orig_files, [input_spline_trace]] )

    print "*"
    print all_ipod_vs_inp_vals
    v6_mat_repmin, v6_mat_repmax = get_replevel_extrema( calc_v6_chipsub_mats, [all_ipod_vs_inp_vals, all_chip_vs_inp_vals])

    #v6_mat_repmin = calc_v6_chipsub_mats(ipod_vs_inp_replo, chip_vs_inp_rephi)
    #v6_mat_repmax = calc_v6_chipsub_mats(ipod_vs_inp_rephi, chip_vs_inp_replo)

    ## write the output info on rzscores
    ci_lo_v6, ci_hi_v6 = get_ci_vecs(v6_mat_boot)

    numpy.save("%s_v6rz_chipsub_orig.npy" % out_prefix, v6_mat_orig)
    numpy.save("%s_v6rz_chipsub_boot.npy" % out_prefix, v6_mat_boot)
    numpy.save("%s_v6rz_chipsub_repmin.npy" % out_prefix, v6_mat_repmin)
    numpy.save("%s_v6rz_chipsub_repmax.npy" % out_prefix, v6_mat_repmax)

    ipod_utils.write_grfile(G_LOCS, v6_mat_orig, "%s_v6rz_chipsub.gr" % out_prefix)
    ipod_utils.write_grfile(G_LOCS, ci_lo_v6, "%s_v6rz_chipsub_95cilo.gr" % out_prefix)
    ipod_utils.write_grfile(G_LOCS, ci_hi_v6, "%s_v6rz_chipsub_95cihi.gr" % out_prefix)
    ipod_utils.write_grfile(G_LOCS, v6_mat_repmin, "%s_v6rz_chipsub_repmin.gr" % out_prefix)
    ipod_utils.write_grfile(G_LOCS, v6_mat_repmax, "%s_v6rz_chipsub_repmax.gr" % out_prefix)

    # also save the same set of gr files for the p-value transformed data
    ipod_utils.write_grfile(G_LOCS, calc_signed_log10p(v6_mat_orig), "%s_v6rzlog10p_chipsub.gr" % out_prefix)
    ipod_utils.write_grfile(G_LOCS, calc_signed_log10p(ci_lo_v6), "%s_v6rzlog10p_chipsub_95cilo.gr" % out_prefix)
    ipod_utils.write_grfile(G_LOCS, calc_signed_log10p(ci_hi_v6), "%s_v6rzlog10p_chipsub_95cihi.gr" % out_prefix)
    ipod_utils.write_grfile(G_LOCS, calc_signed_log10p(v6_mat_repmin), "%s_v6rzlog10p_chipsub_repmin.gr" % out_prefix)
    ipod_utils.write_grfile(G_LOCS, calc_signed_log10p(v6_mat_repmax), "%s_v6rzlog10p_chipsub_repmax.gr" % out_prefix)


