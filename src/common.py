#!/usr/bin/env python

import math
import os
from collections import defaultdict


WORK_DIR = "/Users/agoncear/projects/Enolase/data/"

############### WEIGHTING FUNCTIONS DEPENDING ON THE NUMBERS OF PROTEINS IN EACH COMBINATION ##############
def weight_nprot(N_proteins, combination):
    """ Linear """
    # assuming N_proteins is always > 0
    return N_proteins


def weight_log_nprot(N_proteins, combination):
    """ Log scaling """
    # assuming N_proteins is always > 0
    return math.log(N_proteins)


def weight_zero_one(N_proteins, combination):
    """ Stepwise function 0,1 """
    if N_proteins > len(combination) > 1:
        return 1
    return 0


def weight_one(N_proteins, combination):
    """ Constant 1 """
    return 1


########## FREQUENCIES #################
def calculate_frequencies(fname, weight=weight_one):
    """
    fname: table2.csv
    creates matrices for profile_freq, profile_pair_freq and list of initial clusters
    + def calc_MI
    weight: lambda function
    """

    print("Calculating frequencies for "+fname)
    profile_freq = defaultdict(float)  # key p
    profile_pair = defaultdict(float)  # key (p,q)
    profile_pair_freq = defaultdict(lambda: defaultdict(float))  # [p][q]

    combinations = []  # not used
    labels = []

    label_factors = ['NO']
    with open(fname) as f:
        for line in f:
            # print line
            line = line.strip()
            if line.startswith("Prof") or line.startswith("Total"):
                continue
            # 64,82;92;mandelate racemase: 76_muconate cycloisomerase: 1;10(5)
            comb_line, N_proteins, groupwise, pred = line.split(";", 3)
            label = 0
            groups = {}
            # mandelate racemase: 76_muconate cycloisomerase: 1
            # print groupwise
            if len(groupwise) > 0:
                for group in groupwise.split("_"):
                    # print group
                    g, g_n = group.split(":")
                    g = g.strip()
                    g_n = int(g_n.strip())
                    groups[g] = g_n
                # Only consider the first group as label
                g = groups.keys()[0]
                if g not in label_factors:
                    label_factors.append(g)
                label = label_factors.index(g)
            # 10(5)
            N_in_SF, N_not_in_SF = map(int, pred[:-1].split("(", 1))

            # 92
            N_proteins = int(N_proteins)
            # combination of profiles: 64,82
            comb = set(map(int, comb_line.split(",")))
            if len(comb) < 2:
                continue

            for profile in comb:
                profile_freq[profile] += weight(N_proteins, comb)

            for i, p in enumerate(comb):
                for j, q in enumerate(comb):
                    profile_pair[(p, q)] += weight(N_proteins, comb)
            combinations.append(comb)
            labels.append(label)

    # print label_factors
    P = profile_freq.keys()
    N = len(P)
    S = sum(profile_freq.values()) + N
    profile_freq = {p: (f + 1.0) / S for p, f in profile_freq.iteritems()}

    # normalize pair frequencies
    S = 0.0
    for i, p in enumerate(P):
        for j, q in enumerate(P):
            S += 1.0 + profile_pair[(p, q)]

    for i, p in enumerate(P):
        for j, q in enumerate(P):
            profile_pair_freq[p][q] = (profile_pair[(p, q)] + 1.0) / S  # pseudocount 1.0

    return profile_freq, profile_pair_freq, combinations, labels, label_factors


################# IO PROFILE FREQ ######################
def write_profile_freq(freq, fname):
    print "creating file with profile frequencies (profile_freq.tab): ", fname
    with open(fname, "w") as o:
        for p, f in freq.iteritems():
            o.write("{}\t{}\n".format(p, f))


def get_profile_freq(fname):
    profile_freq = {}
    with open(fname, "r") as f:
        for line in f:
            p, value = line.strip().split(":")
            p = int(p)
            profile_freq[p] = float(value)
    return profile_freq


def write_profile_pair_freq(freq, fname):
    print "creating file with profile pair frequencies (profile_pair_freq.tab)", fname
    with open(fname, "w") as o:
        for p, v in freq.iteritems():
            for q, f in v.iteritems():
                o.write("{}\t{}\t{}\n".format(p, q, f))


def get_profile_pair_freq(fname):
    profile_pair_freq = defaultdict(lambda: defaultdict(float))
    with open(fname, "r") as f:
        for line in f:
            p, q, value = line.strip().split(":", 2)
            p = int(p)
            q = int(q)
            value = float(value)
            profile_pair_freq[p][q] = value
    return profile_pair_freq


################## CYTOSCAPE ####################
def write_cytoscape_edges():
    """
    creates edge file for cytoscape

    input:
    1/MI_matrix.csv - MI matrix of step 1 (before combining two clusters)

    output:
    1/cyt_edge_file.csv - edge file for step 1
    """

    fname = WORK_DIR + "1/MI_matrix.csv"
    sum_mi = 0
    number_of_lines = 0
    with open(fname) as f:
        for line in f:
            number_of_lines += 1
            line = line.strip()
            cl1, cl2, mi = line.split(":", 2)
            mi = float(mi)
            sum_mi += mi
    average_mi = sum_mi/float(number_of_lines)
    with open(fname) as f, open(WORK_DIR + "1/cyt_edge_file.csv", "w") as output_file:
        output_file.write("cluster1\tIC_IC\tcluster2\tMI-score\n")
        for line in f:
            line = line.strip()
            cl1, cl2, mi = line.split(":", 2)
            if float(mi) > average_mi:
                output_file.write("%s\tCL_CL\t%s\t%s\n" % (cl1,cl2,mi))

