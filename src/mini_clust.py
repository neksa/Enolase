#!/usr/bin/env python
"""
1. Calculate frequencies
input:
 * table2.csv - occuring profile combinations in profile-protein_sequence hits

output:
 * profile frequencies (F)
 * profile pair frequencies (G)
 * combinations of profiles (c)
 

2. Default clustering procedure
input:
  * profile frequencies (F)
  * profile pair frequencies (G)

output:
    * n/MI_matrix.csv - MI matrix of step n (before combining two clusters)
    * n/clusters.csv - current set of clusters of step n (before combining two clusters)
    * tree.csv - each line in this file corresponds to two clusters that are going to be combined and MI between them 

3. Fuzzy clustering procedure
input:
 * profile freqs


"""

# import sys
# import math
import random
import operator
import numpy as np

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
# from sklearn.decomposition import PCA

# from matplotlib import pyplot as plt
# from __future__ import print_function

import matplotlib.pyplot as plt
import matplotlib.cm as cm
# from matplotlib.collections import LineCollection
# from mpl_toolkits.mplot3d import Axes3D

from common import *
from clustering import *
from distances import *


def silh(Z, X, clusters):
    range_n_clusters = range(2, 20)

    for n_clusters in range_n_clusters:
        # Create a subplot with 1 row and 2 columns
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(18, 7)

        # The 1st subplot is the silhouette plot
        # The silhouette coefficient can range from -1, 1 but in this example all
        # lie within [-0.1, 1]
        ax1.set_xlim([-0.1, 1])
        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax1.set_ylim([0, len(Z) + (n_clusters + 1) * 10])

        # Initialize the clusterer with n_clusters value and a random generator
        # seed of 10 for reproducibility.
        clusterer = KMeans(n_clusters=n_clusters, random_state=10)
        cluster_labels = clusterer.fit_predict(Z)

        # The silhouette_score gives the average value for all the samples.
        # This gives a perspective into the density and separation of the formed
        # clusters
        silhouette_avg = silhouette_score(Z, cluster_labels)
        print("For n_clusters =", n_clusters,
              "The average silhouette_score is :", silhouette_avg)

        # Compute the silhouette scores for each sample
        sample_silhouette_values = silhouette_samples(Z, cluster_labels)

        y_lower = 10
        for i in range(n_clusters):
            # Aggregate the silhouette scores for samples belonging to
            # cluster i, and sort them
            ith_cluster_silhouette_values = \
                sample_silhouette_values[cluster_labels == i]

            ith_cluster_silhouette_values.sort()

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            color = cm.spectral(float(i) / n_clusters)
            ax1.fill_betweenx(np.arange(y_lower, y_upper),
                              0, ith_cluster_silhouette_values,
                              facecolor=color, edgecolor=color, alpha=0.7)

            # Label the silhouette plots with their cluster numbers at the middle
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples

        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")

        # The vertical line for average silhoutte score of all the values
        ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

        # 2nd Plot showing the actual clusters formed
        colors = cm.spectral(cluster_labels.astype(float) / n_clusters)
        ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30*4, lw=0, alpha=0.7,
                    c=colors)

        # Labeling the clusters
        # centers = clusterer.cluster_centers_
        # # Draw white circles at cluster centers
        # ax2.scatter(centers[:, 0], centers[:, 1],
        #             marker='o', c="white", alpha=1, s=200)

        # for i, c in enumerate(centers):
        #     ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1, s=50)

        ax2.set_title("The visualization of the clustered data.")
        ax2.set_xlabel("Feature space for the 1st feature")
        ax2.set_ylabel("Feature space for the 2nd feature")

        plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                      "with n_clusters = %d" % n_clusters),
                     fontsize=14, fontweight='bold')

        # plt.show()
        fig.savefig("plots/silhouette_{}.png".format(n_clusters), dpi=150)
        plt.clf()


def mini_clustering():
    # hits_fname = WORK_DIR + "table2.csv"
    hits_fname = "/Users/agoncear/projects/Enolase/results/" + "table2.csv"
    F, G, C, L, LF = calculate_frequencies(hits_fname, weight=weight_log_nprot)  # weight_nprot
    superset = reduce(lambda x, y: x | y, C)
    superlist = sorted(list(superset))
    N = len(C)
    D = len(superlist)
    A = np.zeros((N, D))
    for i, c in enumerate(C):
        for j in [superlist.index(p) for p in c]:
            # print i, j, 1
            A[i, j] = 1
    # Remove combinations with only 1 or two profiles
    remlist = set()
    for i, c in enumerate(C):
        if len(c) < 3:  # 3+
            remlist.add(i)
    # Also exclude non-annotated combinations
    for i, l in enumerate(L):
        if l == 0:
            remlist.add(i)

    C = [C[i] for i in range(len(C)) if i not in remlist]
    L = [L[i] for i in range(len(L)) if i not in remlist]
    N = len(C)

    # CALC MATRIX
    seed = np.random.RandomState(seed=3)
    similarity_function = dist00  # dist2  # dist00

    M = np.zeros((N, N))
    for i, c1 in enumerate(C):
        for j, c2 in enumerate(C):
            M[i, j] = similarity_function(c1, c2, F, G)
            # print ",".join(map(str, c1)), ",".join(map(str, c2)), M[i, j]
    print "Shape", M.shape
    print "Min/max", np.min(M), np.max(M)

    with open("plots/initial_similarity_graph.tab", 'w') as o:
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if i < j:
                    ci = "_".join(map(str, sorted(list(C[i]))))
                    cj = "_".join(map(str, sorted(list(C[j]))))
                    o.write("p{}\tp{}\t{}\t{}\t{}\n".format(i, j, ci, cj, M[i, j]))

    C_rand = []
    profiles = F.keys()
    for c in C:
        k = len(c)
        c_rand = set()
        for i in range(k):
            c_rand.add(random.choice(profiles))
        C_rand.append(c_rand)
    M_rand = np.zeros((N, N))
    for i, c1 in enumerate(C_rand):
        for j, c2 in enumerate(C_rand):
            M_rand[i, j] = similarity_function(c1, c2, F, G)
    print "Shape (rand)", M_rand.shape
    print "Min/max (rand)", np.min(M_rand), np.max(M_rand)

    with open("plots/initial_similarity_graph_rand.tab", 'w') as o:
        for i in range(M_rand.shape[0]):
            for j in range(M_rand.shape[1]):
                if i < j:
                    ci = "_".join(map(str, sorted(list(C_rand[i]))))
                    cj = "_".join(map(str, sorted(list(C_rand[j]))))
                    o.write("p{}\tp{}\t{}\t{}\t{}\n".format(i, j, ci, cj, M_rand[i, j]))

    mds = manifold.MDS(n_components=2, metric=True, max_iter=3000, eps=1e-9, dissimilarity="precomputed", random_state=seed)
    X = mds.fit_transform(M)
    print "2D stress", mds.stress_

    fig, ax = plt.subplots()
    ax.scatter(X[:, 0], X[:, 1], c=L)
    # ax.scatter(delta1[:-1], delta1[1:], c=close, s=volume, alpha=0.5)
    ax.set_xlabel(r'Component $\Delta_1$', fontsize=20)
    ax.set_ylabel(r'Component $\Delta_2$', fontsize=20)
    ax.set_title('Reduced space, original labels')
    ax.grid(True)
    # fig.tight_layout()
    # fig1 = plt.gcf()
    fig.savefig("plots/original_labels_in_2D.png", dpi=300)
    plt.clf()

    # model = manifold.TSNE(n_components=3, random_state=seed)
    # X = model.fit_transform(M)

    mds = manifold.MDS(n_components=100, metric=True, max_iter=300, eps=1e-9, dissimilarity="precomputed", random_state=seed)
    Z = mds.fit_transform(M)
    print "Z stress", mds.stress_

    # CLUSTERING KMEANS
    n_clusters = 4
    est = KMeans(n_clusters=n_clusters)
    # est.fit(X)
    est.fit(Z)
    labels = est.labels_
    fig, ax = plt.subplots()
    ax.scatter(X[:, 0], X[:, 1], c=labels)
    # ax.scatter(delta1[:-1], delta1[1:], c=close, s=volume, alpha=0.5)
    ax.set_xlabel(r'Component $\Delta_1$', fontsize=20)
    ax.set_ylabel(r'Component $\Delta_2$', fontsize=20)
    ax.set_title('Kmeans k={} labels in RD100 space'.format(n_clusters))
    ax.grid(True)
    # fig.tight_layout()
    # fig1 = plt.gcf()
    fig.savefig("plots/mds_100_kmeans_{}.png".format(n_clusters), dpi=300)
    plt.clf()

    silh(Z, X, labels)

    # LABELS
    lab_prot = defaultdict(lambda: defaultdict(int))
    # lab_orig_lab = defaultdict(set)
    lab_orig_lab = defaultdict(lambda: defaultdict(int))
    for i in range(M.shape[0]):
        l0 = L[i]
        l = labels[i]
        for c in C[i]:
            lab_prot[l][c] += 1
        # lab_orig_lab[l].add(l0)
        lab_orig_lab[l][l0] += 1
    for l in lab_orig_lab.keys():
        print "Cluster #{} (size N={})".format(l+1, sum(lab_orig_lab[l].values()))
        print "   original labels"
        # , lab_orig_lab[l]
        for k, v in lab_orig_lab[l].iteritems():
            print "             ", LF[k], round(float(v)/sum(lab_orig_lab[l].values()), 2)
        # print "   original labels", lab_orig_lab[l]
        print "   profiles"
        # , lab_prot[l]
        # for k, v in lab_prot[l].iteritems():
        for k, v in sorted(lab_prot[l].items(), key=operator.itemgetter(0)):
            print "             ", k, round(float(v)/sum(lab_orig_lab[l].values()), 2)  # number of labels is the total number of proteins

    # FDR

    # Distance distributions in transformed space
    MM = np.zeros((N, N))
    with open("plots/kmeans_similarity_graph_2D.tab", 'w') as o:
        for i in range(MM.shape[0]):
            for j in range(MM.shape[1]):
                if i < j:
                    ci = "_".join(map(str, sorted(list(C[i]))))
                    cj = "_".join(map(str, sorted(list(C[j]))))
                    MM[i, j] = euclidean_distances(X[i, :], X[j, :])
                    o.write("p{}\tp{}\t{}\t{}\t{}\n".format(i, j, ci, cj, MM[i, j]))

    MM = np.zeros((N, N))
    with open("plots/kmeans_similarity_graph.tab", 'w') as o:
        for i in range(MM.shape[0]):
            for j in range(MM.shape[1]):
                if i < j:
                    ci = "_".join(map(str, sorted(list(C[i]))))
                    cj = "_".join(map(str, sorted(list(C[j]))))
                    MM[i, j] = euclidean_distances(Z[i, :], Z[j, :])
                    o.write("p{}\tp{}\t{}\t{}\t{}\n".format(i, j, ci, cj, MM[i, j]))

    MM_rand = np.zeros((N, N))
    for i, c1 in enumerate(C_rand):
        for j, c2 in enumerate(C_rand):
            MM_rand[i, j] = euclidean_distances(Z[i, :], Z[j, :])
    print "kmeans: Shape (rand)", MM_rand.shape
    print "kmeans: Min/max (rand)", np.min(MM_rand), np.max(MM_rand)

    with open("plots/kmeans_similarity_graph_rand.tab", 'w') as o:
        for i in range(MM_rand.shape[0]):
            for j in range(MM_rand.shape[1]):
                if i < j:
                    ci = "_".join(map(str, sorted(list(C_rand[i]))))
                    cj = "_".join(map(str, sorted(list(C_rand[j]))))
                    o.write("p{}\tp{}\t{}\t{}\t{}\n".format(i, j, ci, cj, MM_rand[i, j]))

"""
library(data.table)
f <- fread("plots/initial_similarity_graph.tab")
g <- fread("plots/initial_similarity_graph_rand.tab")
h <- fread("plots/kmeans_similarity_graph.tab")
x <- fread("plots/kmeans_similarity_graph_2D.tab")

df <- density(f$V5)
dg <- density(g$V5)
dh <- density(h$V5)
dx <- density(x$V5)

plot(range(df$x, dg$x, dh$x, dx$x), range(df$y, dg$y, dh$y, dx$y), type="n", xlab="Distance", ylab="Density")
lines(df, col="black", lwd=2)
lines(dg, col="red", lwd=2, lty=2)
lines(dh, col="green", lwd=2)
lines(dx, col="blue", lwd=2)
legend("topleft", c("Initial dissimilarity", "Random dissimilarity", "Euclidean dist in RD100 space", "Euclidean dist in RD2 space"), text.col = c("black", "red", "green", "blue"))
dev.copy2pdf(file="plots/pairwise_distances.pdf", width=8, height=6)
"""


if __name__ == '__main__':
    np.random.seed(5)
    mini_clustering()

