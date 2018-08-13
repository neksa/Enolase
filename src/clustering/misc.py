

    """
    import skfuzzy as fuzz
    ks = []
    fpcs = []

    from sklearn.cluster import KMeans

    for k in range(2, 200):
        kms = KMeans(k)
        kms.fit(A.T)
        print kms.get_params()
        ks.append(k)
        fpcs.append(kms.score)
        
        # cntr, U, U0, d, Jm, p, fpc = fuzz.cluster.cmeans(A.T, k, 2.0, error=0.05, maxiter=1000)  # , U_init=None)
        # print cntr
        # print cntr.shape
        # print Jm
        # print k, fpc
        # ks.append(k)
        # fpcs.append(fpc)
        # print U
        # print cntr, U, U0, d, Jm, p, fpc

    cntr, U, U0, d, Jm, p, fpc = fuzz.cluster.cmeans(A.T, 2, 2.0, error=0.005, maxiter=1000, U_init=None)
    print cntr
    """


    """
    # import matplotlib.pyplot as plt
    import pylab as pl
    pl.clf()
    pl.plot(ks, fpcs, 'k--', color="blue")
    pl.xlabel('K clusters')
    pl.ylabel('Fuzzy partition coefficient, fpc')
    pl.title("")
    # pl.legend(loc="lower right")
    pl.savefig("k_vs_fpc.png")
    """

    """
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    pca.fit(A)
    print(pca.explained_variance_ratio_)
    X = pca.transform(A)
    # print X

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig, ax = plt.subplots()
    # print len(X[:,0]), len(L)
    ax.scatter(X[:, 0], X[:, 1], c=L)
    # ax.scatter(delta1[:-1], delta1[1:], c=close, s=volume, alpha=0.5)
    ax.set_xlabel(r'$\Delta_1$', fontsize=20)
    ax.set_ylabel(r'$\Delta_2$', fontsize=20)
    ax.set_title('Proteins clustered by profile combinations')
    ax.grid(True)
    # fig.tight_layout()
    # fig1 = plt.gcf()
    fig.savefig("components.png", dpi=300)

    fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    ax = Axes3D(fig)
    pca3 = PCA(n_components=3)
    pca3.fit(A)
    print(pca3.explained_variance_ratio_)
    Z = pca3.transform(A)
    ax.scatter(Z[:, 0], Z[:, 1], Z[:, 2], c=L)
    ax.grid(True)
    fig.savefig("components_3D.png", dpi=300)
    """