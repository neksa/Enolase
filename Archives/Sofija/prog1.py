#!/usr/bin/env python
"""
clustering procedure

input:
no input, but program prog4.py shoul be run first

output:
n/MI_matrix.csv - MI matrix of step n (before combining two clusters)
n/clusters.csv - current set of clusters of step n (before combining two clusters)
tree.csv - each line in this file corresponds to two clusters that are going to be combined and MI between them 
"""
import os
import math
import prog4
import common
WORK_DIR = common.WORK_DIR

profile_freq = prog4.get_profile_freq()
profile_pair_freq = prog4.get_profile_pair_freq()
clusters = prog4.get_clusters()


def calc_MI_matrix():
    matrix = {}
    for i, c1 in enumerate(clusters):
        for j, c2 in enumerate(clusters):
            if j >= i: continue
            matrix[(i, j)] = round(prog4.calc_MI_fast(c1, c2, profile_pair_freq, profile_freq), 4)
    return matrix

def create_MI_file_for_this_step():
    if not os.path.exists(WORK_DIR + str(step)):
        os.makedirs(WORK_DIR + str(step))
    nf = open(WORK_DIR + str(step) + "/MI_matrix.csv","w")
    for p,v in C.iteritems():
        #print v
        nf.write("%s:%s:%.5f\n" % (",".join(map(str,sorted(list(map(int,clusters[p[0]]))))), ",".join(map(str,sorted(list(map(int,clusters[p[1]]))))),v))
    nf.close()

def create_file_with_information_about_each_cluster_for_this_step():
    if not os.path.exists(WORK_DIR + str(step)):
        os.makedirs(WORK_DIR + str(step))
    nf = open(WORK_DIR + str(step) + "/clusters.csv","w")
    for c in clusters:
        nf.write(",".join(map(str,sorted(list(map(int,c))))) + '\n')
    nf.close()

tree_file = open(WORK_DIR + "tree.csv","w")
step = 1
C = {}  # Mutual information scores for pairs of clusters
while (len(clusters) > 1):
    print "step " + str(step)
    C = calc_MI_matrix()
    pair, value = max(C.iteritems(), key=lambda x:x[1])
    print clusters[pair[0]],clusters[pair[1]],value
    
    combining_identical_clusters = False
    for a, c1 in enumerate(clusters):
        for b, c2 in enumerate(clusters):
            if a>b and c1.issubset(c2) and c2.issubset(c1) and a != b:
                pair = a,b
                combining_identical_clusters = True
                break;
    i, j = pair
    if not combining_identical_clusters:
        tree_file.write("%s:%s:%.5f\n" % (",".join(map(str,sorted(list(map(int,clusters[i]))))), ",".join(map(str,sorted(list(map(int,clusters[j]))))), value))
        create_MI_file_for_this_step()
        create_file_with_information_about_each_cluster_for_this_step()
        step += 1
    clusters.append(clusters[i] | clusters[j])
    clusters = [cl for k, cl in enumerate(clusters) if k not in (i, j)]
tree_file.close()
C = calc_MI_matrix()
create_MI_file_for_this_step()
create_file_with_information_about_each_cluster_for_this_step()
