import math

"""
Comparing profile combinations
"""

def calc_all_freqs():
    hits_fname = WORK_DIR + "table2.csv"
    for weight, weight_fname in zip((weight_one, weight_zero_one, weight_nprot, weight_log_nprot), ('one', 'zero_one', 'nprot', 'log_nprot')):
        F, G, C, L = calculate_frequencies(hits_fname, weight=weight)
        F_file = WORK_DIR + "F_{}.tab".format(weight_fname)
        G_file = WORK_DIR + "G_{}.tab".format(weight_fname)
        write_profile_freq(F, F_file)
        write_profile_pair_freq(G, G_file)


def calc_agglomerative_clustering():
    hits_fname = WORK_DIR + "table2.csv"
    F, G, C, L = calculate_frequencies(hits_fname, weight=weight_log_nprot)
    # F = get_profile_freq(WORK_DIR + 'F_log_nprot.tab')
    # G = get_profile_pair_freq(WORK_DIR + 'G_log_nprot.tab')
    agglomerative_clustering(F, G, C)


# pointwise MI normalized by the length of profiles
def dist3(c1, c2, F, G):
    MI = 0.0
    for p in c1:
        for q in c2:
            if G[p][q] > 0:
                NF = - math.log(G[p][q], 2)
                MI += math.log(G[p][q] / (F[p] * F[q]), 2) / NF

    MI = MI / (len(c1) * len(c2))
    return MI


# normalized pointwise MI
def dist2(c1, c2, F, G):
    MI = 0.0
    for p in c1:
        for q in c2:
            if G[p][q] > 0:
                NF = - math.log(G[p][q], 2)
                MI += math.log(G[p][q] / (F[p] * F[q]), 2) / NF
    return MI


# pure pointwise mutual information
def dist1(c1, c2, F, G):
    MI = 0.0
    for p in c1:
        for q in c2:
            if G[p][q] > 0.0:
                MI += math.log(G[p][q] / (F[p] * F[q]), 2)
    # return 200.0 - MI
    return MI


# joint entropy
def H(c1, c2, F, G):
    h = 0.0
    for p in c1:
        for q in c2:
            # if p == q:
            #     h -= F[p] * math.log(F[p], 2)
            # else:
            if G[p][q] > 0.0:
                h -= G[p][q] * math.log(G[p][q], 2)
    return h


# pure mutual information
def MI(c1, c2, F, G):
    mi = 0.0
    for p in c1:
        for q in c2:
            # if p == q:
            #     mi += F[p] * math.log(F[p], 2)
            # else:
            if G[p][q] > 0.0:
                # print p, q, G[p][q], F[p], F[q], G[p][q] * math.log(G[p][q] / (F[p] * F[q]), 2)
                mi += G[p][q] * math.log(G[p][q] / (F[p] * F[q]), 2)
    # print c1, c2, mi
    return mi


# variation of information
# Compares clusterings, is a metric
def VI(c1, c2, F, G):
    mi = MI(c1, c2, F, G)
    h1 = 0
    for p in c1:
        h1 -= F[p] * math.log(F[p], 2)
    h2 = 0
    for q in c2:
        h2 -= F[q] * math.log(F[q], 2)

    vi = h1 + h2 - 2*mi
    return vi


# normalized distance based on mutual information
def dist00(c1, c2, F, G):
    h = H(c1, c2, F, G)
    mi = MI(c1, c2, F, G)
    # print "dist00", c1, c2, h, mi
    if h == 0.0:
        return 0.0
    d = h - mi
    return d / h
