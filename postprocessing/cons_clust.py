#!/usr/bin/python

# perform consensus clustering, based on the philosophy of iclust (doi: 10.1073/pnas.0507432102)

import scipy.cluster.vq as vc
import scipy.cluster.hierarchy as hc
import scipy.spatial as sp
import sklearn.cluster
import numpy

def cons_clust( datmat, ntrials, nclust, max_clust, plotfile=None):
    # perform consensus clustering on an input data matrix
    # arguments:
    #   datmat: the input data matrix, with observations to be clustered as rows
    #   ntrials: number of kmeans runs to do AT EACH VALUE OF NCLUST
    #   nclust: either integer or list of integers. number of clusters to use in the k-means step
    #   max_clust: maximum number of final clusters to be returned
    #   plotfile: an optional plot file for the co-clustering matrix
    # returns:
    #   a list of final cluster ids for each observation

    n_vars = datmat.shape[0]
    together_count = numpy.zeros( (n_vars, n_vars), dtype='int' )

    if type(nclust) == int:
        nclust_list = [nclust]
    else:
        nclust_list = nclust

    # now do the specified number of clustering trials, recording how often different observations go together
    for n_c in nclust_list:
        for iter_num in range(ntrials):
            k_fit = sklearn.cluster.KMeans(n_clusters=n_c).fit(datmat)
            assignments = k_fit.labels_

            for i in range(n_vars):
                for j in range(i,n_vars):
                    if assignments[i] == assignments[j]:
                        together_count[i,j] += 1
                        if i != j:
                            together_count[j,i] += 1

    together_frac = together_count / float(together_count[0,0])
    #DEBUG print "yyy"
    #DEBUG print float(together_count[0,0])
    #DEBUG print together_frac
    #DEBUG print "zzz"


    # now, we can do linkage clustering based on this fraction of the time that observations go together
    cond_dist_mat = sp.distance.squareform( 1.0 - together_frac )
    clust_obj = hc.complete(cond_dist_mat)
    print cond_dist_mat
    cluster_ids = hc.fcluster( clust_obj, max_clust, criterion='maxclust' )
    if plotfile is not None:
        import seaborn
        import matplotlib.pyplot as plt
        plt.figure()
        seaborn.clustermap( together_frac, row_linkage = clust_obj, col_linkage=clust_obj, cmap='Blues', vmin=0.0, vmax=1.0, square=True)
        plt.savefig(plotfile)
    #alternative -- use a threshold: cluster_ids = hc.fcluster( clust_obj, 0.4, criterion='inconsistent' )
    return cluster_ids



