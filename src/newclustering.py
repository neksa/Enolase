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

import sys
import numpy as np

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D

from common import *
from clustering import *


def new_clustering():
    hits_fname = WORK_DIR + "table2.csv"
    F, G, C, L = calculate_frequencies(hits_fname, weight=weight_log_nprot)
    # F, G, C, L = calculate_frequencies(hits_fname, weight=weight_nprot)
    superset = reduce(lambda x, y: x | y, C)
    superlist = sorted(list(superset))
    N = len(C)
    D = len(superlist)
    A = np.zeros((N, D))
    # print A
    # print superlist.index(0)
    for i, c in enumerate(C):
        for j in [superlist.index(p) for p in c]:
            # print i, j, 1
            A[i, j] = 1

    # print A
    # print A.shape
    # print L

    # heat kernel: distance -> similarity
    # in case of signed distance matrix, etc
    # similarity = np.exp(-beta * distance / distance.std())
    # Gaussian kernel
    #  np.exp(-gamma * d(X,X) ** 2)

    # .SpectralClustering

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
    # similarity_function = VI  # dist2  # dist00

    M = np.zeros((N, N))
    for i, c1 in enumerate(C):
        for j, c2 in enumerate(C):
            M[i, j] = similarity_function(c1, c2, F, G)
            print ",".join(map(str, c1)), ",".join(map(str, c2)), M[i, j]
    print "Shape", M.shape
    print "Min/max", np.min(M), np.max(M)

    # Inverse:
    # M = M - np.min(M)
    # M = 0 - M
    # print "New min/max", np.min(M), np.max(M)

    with open("similarity_graph.tab", 'w') as o:
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if i < j:
                    ci = "_".join(map(str, sorted(list(C[i]))))
                    cj = "_".join(map(str, sorted(list(C[j]))))
                    o.write("p{}\tp{}\t{}\t{}\t{}\n".format(i, j, ci, cj, M[i, j]))

    import random
    # random.randint(3, 8)
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

    with open("similarity_graph_rand.tab", 'w') as o:
        for i in range(M_rand.shape[0]):
            for j in range(M_rand.shape[1]):
                if i < j:
                    ci = "_".join(map(str, sorted(list(C_rand[i]))))
                    cj = "_".join(map(str, sorted(list(C_rand[j]))))
                    o.write("p{}\tp{}\t{}\t{}\t{}\n".format(i, j, ci, cj, M_rand[i, j]))

    # sys.exit(0)



    # x = []
    # y = []
    # for N in range(1, 21):
    #     mds = manifold.MDS(n_components=N, metric=True, max_iter=3000, eps=1e-9, random_state=seed,
    #                        dissimilarity="precomputed")  #, n_jobs=1, n_init=1)
    #     X = mds.fit_transform(M)  # distM  #, init=pos)
    #     print X.shape, mds.stress_
    #     x.append(X.shape[1])
    #     y.append(mds.stress_)
    # line, = plt.plot(x, y, '-', linewidth=2)
    # # dashes = [10, 5, 100, 5] # 10 points on, 5 off, 100 on, 5 off
    # # line.set_dashes(dashes)
    # plt.xlabel("N dimensions")
    # plt.ylabel("Stress")
    # plt.savefig("scree.png", dpi=300)


    # model = manifold.TSNE(n_components=2, random_state=0)
    # X = model.fit_transform(M)

    mds = manifold.MDS(n_components=2, metric=True, max_iter=3000, eps=1e-9, dissimilarity="precomputed", random_state=seed)
    X = mds.fit_transform(M)
    print "2D stress", mds.stress_

    fig, ax = plt.subplots()
    ax.scatter(X[:, 0], X[:, 1], c=L)
    # ax.scatter(delta1[:-1], delta1[1:], c=close, s=volume, alpha=0.5)
    ax.set_xlabel(r'$\Delta_1$', fontsize=20)
    ax.set_ylabel(r'$\Delta_2$', fontsize=20)
    ax.set_title('Proteins clustered by profile combinations')
    ax.grid(True)
    # fig.tight_layout()
    # fig1 = plt.gcf()
    fig.savefig("mds.png", dpi=300)
    plt.clf()

    # model = manifold.TSNE(n_components=3, random_state=seed)
    # X = model.fit_transform(M)
 
    mds = manifold.MDS(n_components=3, metric=True, max_iter=3000, eps=1e-9, dissimilarity="precomputed", random_state=seed)
    Z = mds.fit_transform(M)
    print "3D stress", mds.stress_

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(Z[:, 0], Z[:, 1], Z[:, 2], c=L)
    ax.grid(True)
    fig.savefig("mds3D.png", dpi=300)
    plt.clf()


    mds = manifold.MDS(n_components=100, metric=True, max_iter=300, eps=1e-9, dissimilarity="precomputed", random_state=seed)
    Z = mds.fit_transform(M)
    print "Z stress", mds.stress_

    # CLUSTERING KMEANS
    from sklearn.cluster import KMeans
    est = KMeans(n_clusters=6)
    # est.fit(X)
    est.fit(Z)
    labels = est.labels_
    fig, ax = plt.subplots()
    ax.scatter(X[:, 0], X[:, 1], c=labels)
    # ax.scatter(delta1[:-1], delta1[1:], c=close, s=volume, alpha=0.5)
    ax.set_xlabel(r'$\Delta_1$', fontsize=20)
    ax.set_ylabel(r'$\Delta_2$', fontsize=20)
    ax.set_title('Proteins clustered by profile combinations')
    ax.grid(True)
    # fig.tight_layout()
    # fig1 = plt.gcf()
    fig.savefig("mds_clust_kmeans.png", dpi=300)
    plt.clf()


    # CLUSTERING AGGLOMERATIVE
    from sklearn.cluster import AgglomerativeClustering
    aggl = AgglomerativeClustering(n_clusters=6, affinity='precomputed', linkage='average')
    labels = aggl.fit_predict(M)

    # aggl = AgglomerativeClustering(n_clusters=6)
    # labels = aggl.fit_predict(Z)

    # labels = est.labels_
    fig, ax = plt.subplots()
    ax.scatter(X[:, 0], X[:, 1], c=labels)
    # ax.scatter(delta1[:-1], delta1[1:], c=close, s=volume, alpha=0.5)
    ax.set_xlabel(r'$\Delta_1$', fontsize=20)
    ax.set_ylabel(r'$\Delta_2$', fontsize=20)
    ax.set_title('Proteins clustered by profile combinations')
    ax.grid(True)
    # fig.tight_layout()
    # fig1 = plt.gcf()
    fig.savefig("mds_clust_aggl_.png", dpi=300)
    plt.clf()

    from scipy.spatial.distance import cdist, euclidean
    def dunn(points, labels):
        clusters = list(set(labels))
        intra_dist = []
        for c in clusters:
            d = []
            for i, ci in enumerate(labels):
                if ci == c:
                    x = points[i, :]
                    for j, cj in enumerate(labels):
                        if cj == c:
                            if j > i:
                                y = points[j, :]
                                d.append(euclidean(x, y))
            a = np.mean(d)
            if len(d) == 0:
                a = 0.0
            # print "intra", c, a
            intra_dist.append(a)
        max_a = np.max(intra_dist)
        inter_dist = []
        for pi, p in enumerate(clusters):
            for qj, q in enumerate(clusters):
                d = []
                if qj <= pi:
                    continue
                for i, ci in enumerate(labels):
                    if ci == p:
                        x = points[i, :]
                        for j, cj in enumerate(labels):
                            if cj == q:
                                y = points[j, :]
                                d.append(euclidean(x, y))
                b = np.mean(d)
                # print "inter", p, q, b
                inter_dist.append(b)
        min_b = np.min(inter_dist)
        if max_a == 0.0:
            D = 0.0
        else:
            D = min_b / max_a
        return D

    # benchmarking clustering
    x = []
    y = []
    d = []
    for k in range(2, 50):
        est = KMeans(n_clusters=k)
        est.fit(Z)

        labels = est.labels_
        D = dunn(Z, labels)
        # labels = est.labels_
        inertia = est.inertia_
        print inertia, D

        # write plot for k-means iteration
        if k < 12:
            fig, ax = plt.subplots()
            ax.scatter(X[:, 0], X[:, 1], c=labels)
            ax.set_xlabel(r'$\Delta_1$', fontsize=20)
            ax.set_ylabel(r'$\Delta_2$', fontsize=20)
            ax.set_title('K-means with k={}, Dunn={}'.format(k, D))
            ax.grid(True)
            fig.savefig("mds_clust_kmeans_{}.png".format(k), dpi=300)
            plt.clf()

        x.append(k)
        y.append(inertia)
        d.append(D)
    line, = plt.plot(x, y, '-', linewidth=2)
    plt.xlabel("K clusters")
    plt.ylabel("Inertia")
    plt.savefig("clusters_inertia.png", dpi=300)
    plt.clf()

    line, = plt.plot(x, d, '-', linewidth=2)
    plt.xlabel("K clusters")
    plt.ylabel("Dunn index")
    plt.savefig("clusters_dunn.png", dpi=300)
    plt.clf()

    # CLUSTERING FUZZY:
    import skfuzzy as fuzz
    cntr, U, U0, d, Jm, p, fpc = fuzz.cluster.cmeans(Z.T, 6, 2.0, error=1e-9, maxiter=3000)  #, U_init=None)
    print "Fuzzy FPC", fpc
    # Harden Fuzzy cluster membership:
    labels = U.argmax(axis=0)
    print U.shape, labels.shape, labels
    fig, ax = plt.subplots()
    ax.scatter(X[:, 0], X[:, 1], c=labels)
    ax.set_xlabel(r'$\Delta_1$', fontsize=20)
    ax.set_ylabel(r'$\Delta_2$', fontsize=20)
    ax.set_title('Proteins clustered by profile combinations')
    ax.grid(True)
    fig.savefig("mds_clust_fuzzy.png", dpi=300)
    plt.clf()
    # print cntr
    # defuzz_centroid = fuzz.defuzz(x, mfx, 'centroid')  # Same as skfuzzy.centroid
    # PREDICT MEMBERSHIP FOR FUZZY CLUSTERING
    # U, _, _, _, _, fpc2 = fuzz.cluster.cmeans_predict(Z, cntr, 2., error=0.005, maxiter=1000)
    # print "Fuzzy predict FPC", fpc2
    # print U



if __name__ == '__main__':
    # calc_all_freqs()
    # calc_agglomerative_clustering()
    new_clustering()

