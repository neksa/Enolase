
def get_dunn_index(fdist, *clusters):
    """
    Returns the Dunn index for the given selection of nodes.

    J.C. Dunn. Well separated clusters and optimal fuzzy
    partitions. 1974. J.Cybern. 4. 95-104.

    """

    if len(clusters)<2:
        raise ValueError, "At least 2 clusters are required"

    intra_dist = []
    for c in clusters:
        for i in c.get_leaves():
            if i is not None:
                # item intraclsuterdist -> Centroid Diameter
                a = fdist(i.profile, c.profile)*2
                intra_dist.append(a)
    max_a = numpy.max(intra_dist)
    inter_dist = []
    for i, ci in enumerate(clusters):
        for cj in clusters[i+1:]:
            # intracluster dist -> Centroid Linkage
            b = fdist(ci.profile, cj.profile)
            inter_dist.append(b)
    min_b = numpy.min(inter_dist)

    if max_a == 0.0:
        D = 0.0
    else:
        D = min_b / max_a
    return D

