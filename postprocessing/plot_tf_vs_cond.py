#!/usr/bin/python

# given a pickle file containing a database of TFBS occupancies, plot heat maps of the resulting data
# The input pickle file must contain the following columns:
#   condition -- a label describing the biological condition
#   startloc -- the 0-indexed start location of the binding site
#   endloc -- the 0-indexed end location of the binding site 
#   occupancy -- a scalar value giving the occupancy at the site of interest
#   regulator -- a label for the TF under consideration
#   site_ID -- A string giving a more detailed description of the site (for labeling purposes). Often we use TFname_startloc_endloc, although there are no particular requirements

# This pickle file should be specified as the `pkl_file` argument to this script. In addition, at a minimum one must specify the 
#  prefix to use for output, and the minimum number of sites needed for a TF to be included in clustering/plotting. 
# Several other arguments have key impacts on the program's behavior:
#  --conds: indicates the condition names that will be considered in the occupancy plotting and clustering
#  --refcond: indicates the reference condition to use in the difference maps that are provided in output
#  --nclust: target number of clusters to divide data into

# In order to perform clustering, the site-level occupancy data for each TF from the input table will first be aggregated into an
# average occupancy for each condition, and then the profiles of occupancy across conditions for each TF will be clustered.
# The values for all sites for a given TF are usually aggregated using the geometric mean of occupancy scores, but the arithmetic
# mean can be used instead if the --meanavg argument is set

# In addition, the --scalemethod argument indicates how data should be scaled across conditions within each TF to allow appropriate normalization.
#    options are:
#        (the first four options act on site level occupancy, prior to averaging the values within each condition)
#        identity -- no scaling
#        mean -- divide each site-value by the mean site-level occupancy for that TF across conditions, after having subtracted the minimum observed site-level value for that TF
#        max -- divide each site-value by the max site-level occupancy for that TF across conditions, after having subtracted the minimum observed site-level value for that TF
#        iqr -- subtract from each site-level value the median site-level value for that TF, and divide by the iqr for observed sites for that regulator
#        (the last three options act on condition level occupancy, after averaging the values within each condition)
#        maxcond -- divide the condition-level average occupancy for each TF by the value for the highest site level average occupancy for that TF
#        medcond -- divide the value for each condition by the median condition-level value for that TF
#        rankcond -- assign numerical ranks from 1 to n (for n conditions) for the condition-level average values for each TF

# The key output files generated have the following suffixes:
#_unsorted.pdf: heat map of average occupancy under each condition prior to clustering
#_sorted.pdf: heat map as in _unsorted.pdf, but with the values sorted by the clustering results
#_clust_ids.pdf: A pdf showing colors corresponding to the cluster IDs, which can be placed directly alongside the _sorted.pdf image to show the cluster IDs
#_occs_for_plot.csv: Table containing the normalized average occupancy for each TF under each condition, plus the cluster IDs that were calculated
#Also, files with the same suffixes as above but also containing _diffs will be generated; these correspond to the application of the same clustering procedure as noted above, but on the DIFFERENCES between each condition-level average score and the score from a reference condition

    parser.add_argument('--conds',help="conditions to use (default: use all conditions in pkl file)", default=None)
    parser.add_argument('--nclust',help="number of clusters to use (default 10)", default=10, type=int)
    parser.add_argument('--meanavg',help="use arithmetic instead of geometric means", default=False, action='store_true')
    parser.add_argument('--scalemethod',help="method for scaling each TF (options: identity, mean, max, iqr,maxcond,medcond,rankcond)", default="identity")
    parser.add_argument('pkl_file', help='pandas pickle file containing the site-level occupancy data to consider')
    parser.add_argument('outprefix',help="prefix for output files")
    parser.add_argument('mincount',help="Minimum number of binding sites to consider a tf", type=int)
    parser.add_argument('--matfile',help="file to plot the togetherness matrix from clustering", default=None)
    parser.add_argument('--cmap',help="Seaborn color map to use", default="Blues")
    parser.add_argument('--refcond',help="name of the reference condition", default="rdmglu")
    parser.add_argument('--clustreps',help="number of replicates to use for consensus clustering", type=int, default=200)


import matplotlib
matplotlib.use('Agg')

import sys
import argparse
import numpy
import scipy.stats
import scipy.cluster
import matplotlib.pyplot as plt
import seaborn as sns
from bisect import bisect_left
import pandas
import sklearn.cluster
import cons_clust

sns.set()

if __name__ == "__main__":

    # set up commaand line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('--conds',help="conditions to use (default: use all conditions in pkl file)", default=None)
    parser.add_argument('--nclust',help="number of clusters to use (default 10)", default=10, type=int)
    parser.add_argument('--meanavg',help="use arithmetic instead of geometric means", default=False, action='store_true')
    parser.add_argument('--scalemethod',help="method for scaling each TF (options: identity, mean, max, iqr,maxcond,medcond,rankcond)", default="identity")
    parser.add_argument('pkl_file', help='pandas pickle file containing the site-level occupancy data to consider')
    parser.add_argument('outprefix',help="prefix for output files")
    parser.add_argument('mincount',help="Minimum number of binding sites to consider a tf", type=int)
    parser.add_argument('--matfile',help="file to plot the togetherness matrix from clustering", default=None)
    parser.add_argument('--cmap',help="Seaborn color map to use", default="Blues")
    parser.add_argument('--refcond',help="name of the reference condition", default="rdmglu")
    parser.add_argument('--clustreps',help="number of replicates to use for consensus clustering", type=int, default=200)

    args=parser.parse_args()

    tf_db = pandas.read_pickle(args.pkl_file)

    if args.conds is None:
        cond_list = set(tf_db.condition)
    else:
        cond_list = args.conds.split(",")

    tf_db = tf_db.loc[ [x in cond_list for x in tf_db.condition ] ]

    # prune the tf list based on number of sites
    badkeys = []
    for tf in set(tf_db.regulator):
        numsites = len(set(tf_db[tf_db.regulator == tf].site_ID))
        if numsites< args.mincount:
            badkeys += [tf]

    tf_db_pruned = pandas.DataFrame( tf_db.loc[[ not(x in badkeys) for x in tf_db.regulator]] )

    # now do any rescaling by TF that is requested
    if args.scalemethod == "max":
        tf_db_pruned["occupancy"] = numpy.fmax(0.0, tf_db_pruned.occupancy)
        max_by_tf = tf_db_pruned.pivot_table(index='regulator',values='occupancy', aggfunc=max)
        min_by_tf = tf_db_pruned.pivot_table(index='regulator',values='occupancy', aggfunc=min)
        tf_db_pruned["occupancy"] = (tf_db_pruned.occupancy - min_by_tf.loc[tf_db_pruned.regulator].occupancy.values) / (max_by_tf.loc[tf_db_pruned.regulator].occupancy.values - min_by_tf.loc[tf_db_pruned.regulator].occupancy.values)

    if args.scalemethod == "mean":
        avg_by_tf = tf_db_pruned.pivot_table(index='regulator',values='occupancy')
        min_by_tf = tf_db_pruned.pivot_table(index='regulator',values='occupancy', aggfunc=min)
        tf_db_pruned["occupancy"] = (tf_db_pruned.occupancy - min_by_tf.loc[tf_db_pruned.regulator].occupancy.values) / avg_by_tf.loc[tf_db_pruned.regulator].occupancy.values

    if args.scalemethod == "iqr":
        meds_by_tf = tf_db_pruned.pivot_table(index='regulator',values='occupancy', aggfunc=lambda x : numpy.median(x))
        iqr_low = tf_db_pruned.pivot_table(index='regulator',values='occupancy', aggfunc=lambda x : scipy.stats.scoreatpercentile(x, 25))
        iqr_hi = tf_db_pruned.pivot_table(index='regulator',values='occupancy', aggfunc=lambda x : scipy.stats.scoreatpercentile(x, 75))
        iqr = iqr_hi - iqr_low
        tf_db_pruned["occupancy"] = (tf_db_pruned.occupancy - meds_by_tf.loc[tf_db_pruned.regulator].occupancy.values) / iqr.loc[tf_db_pruned.regulator].occupancy.values



    # generate appropriate summary statistics
    if args.meanavg:
        occ_mat = tf_db_pruned.pivot_table(index='regulator',columns='condition', values='occupancy', aggfunc=numpy.median)
    else:
        occ_mat = tf_db_pruned.pivot_table(index='regulator',columns='condition', values='occupancy', aggfunc= (lambda x : scipy.stats.mstats.gmean(numpy.fmax(x,0.01))) )

    if args.scalemethod == 'maxcond':
        print occ_mat
        occ_mat_rowmax = numpy.amax( occ_mat, axis=1)
        occ_mat= occ_mat / numpy.tile(occ_mat_rowmax.values.reshape(-1,1), occ_mat.shape[1])
        print occ_mat_rowmax
        print occ_mat

    if args.scalemethod == 'medcond':
        print occ_mat
        occ_mat_rowmax = numpy.median( occ_mat.values, axis=1)
        occ_mat= occ_mat / numpy.tile(occ_mat_rowmax.reshape(-1,1), occ_mat.shape[1])
        print occ_mat_rowmax
        print occ_mat

    if args.scalemethod == 'rankcond':
        print occ_mat
        occ_mat_ranked = scipy.stats.mstats.rankdata( occ_mat.values, axis=1)
        occ_mat /= occ_mat
        occ_mat *= occ_mat_ranked
        print occ_mat

    #numpy.save('tmpmat.npy', occ_mat)
    f, ax = plt.subplots(figsize=(9, 24))
    sns.heatmap(occ_mat, cmap=args.cmap,robust=True)
    plt.savefig("%s_unsorted.pdf" % args.outprefix)
    print args.nclust
    print min(args.clustreps, occ_mat.shape[0])
    # the commented out line uses standard consensus clustering
    k_lab = cons_clust.cons_clust( occ_mat, 200, range(args.nclust - 2,args.nclust + 2), args.nclust, plotfile=args.matfile)
    #k_lab = cons_clust_gmm.cons_clust( occ_mat, args.clustreps, args.nclust, args.nclust, plotfile=args.matfile)
    print "===="
    print k_lab
    occ_mat["cluster_id"] = k_lab
    occ_mat.to_csv('%s_occs_for_plot_withclust.csv' % args.outprefix)
    occ_mat.sort_values(by="cluster_id", axis=0,inplace=True)
    clust_order = occ_mat["cluster_id"]
    cluster_order_f = pandas.DataFrame( {'cluster_id' : clust_order, 'clust2' : clust_order} )
    occ_mat.drop(labels="cluster_id", axis=1,inplace=True)
    #row_l_mat = numpy.equal.outer(k_lab, k_lab)
    f, ax = plt.subplots(figsize=(9, 40))
    #print "%%%"
    #print row_l_mat
    #print "xxx"
    sns.heatmap(occ_mat, cmap=args.cmap,robust=True,yticklabels=True, vmin=0.0)
    f.set_size_inches(9,40)
    plt.savefig("%s_sorted.pdf" % args.outprefix)

    f, ax = plt.subplots(figsize=(9, 40))
    print clust_order
    sns.heatmap(cluster_order_f, cmap=sns.color_palette("bright"),robust=False,yticklabels=True)
    plt.savefig("%s_clust_ids.pdf" % args.outprefix)

    occ_mat.to_csv('%s_occs_for_plot.csv' % args.outprefix)

    # also try a completely different procedure, where we plot the changes relative to the reference condition rather than the absolute occupancies

    if args.meanavg:
        occ_mat = tf_db_pruned.pivot_table(index='regulator',columns='condition', values='occupancy')
    else:
        occ_mat = tf_db_pruned.pivot_table(index='regulator',columns='condition', values='occupancy', aggfunc= (lambda x : scipy.stats.mstats.gmean(numpy.fmax(x,0.01))) )


    diff_mat = pandas.DataFrame(occ_mat)

    norm_col = numpy.array(occ_mat[args.refcond])
    for c in diff_mat.columns:
        diff_mat.loc[:,c] = occ_mat.loc[:,c] - norm_col

    diff_mat.drop(labels=args.refcond, axis=1, inplace=True)

    f, ax = plt.subplots(figsize=(9, 24))
    sns.heatmap(diff_mat, cmap='RdBu',robust=True, center=0.0, vmin=-1.0, vmax=1.0)
    plt.savefig("%s_diffs_unsorted.pdf" % args.outprefix)
    k_lab = cons_clust.cons_clust( diff_mat, 200, range(args.nclust - 2,args.nclust + 2), args.nclust, plotfile = args.matfile + "_diffmat.pdf")

    diff_mat["cluster_id"] = k_lab
    diff_mat.to_csv('%s_diffs_for_plot_withclust.csv' % args.outprefix)
    diff_mat.sort_values(by="cluster_id", axis=0,inplace=True)
    clust_order = diff_mat["cluster_id"]
    cluster_order_f = pandas.DataFrame( {'cluster_id' : clust_order, 'clust2' : clust_order} )
    diff_mat.drop(labels="cluster_id", axis=1,inplace=True)

    f, ax = plt.subplots(figsize=(9, 40))

    sns.heatmap(diff_mat, cmap='RdBu',robust=True,yticklabels=True, center=0.0, vmin=-2.0, vmax=2.0)
    f.set_size_inches(9,40)
    plt.savefig("%s_diffs_sorted.pdf" % args.outprefix)

    f, ax = plt.subplots(figsize=(9, 40))
    print clust_order
    sns.heatmap(cluster_order_f, cmap=sns.color_palette("bright"),robust=False,yticklabels=True)
    plt.savefig("%s_diffs_clust_ids.pdf" % args.outprefix)

    diff_mat.to_csv('%s_diffs_for_plot.csv' % args.outprefix)
