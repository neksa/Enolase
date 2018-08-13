#!/usr/bin/env python
"""
experimantal clustering procedure

What is different from prog1.py: procedure remove_orphans() which is called in line 78

input:
table2.csv - set of initial clusters

output:
Nexp/MI_matrix.csv - MI matrix of step n (before combining two clusters)
Nexp/clusters.csv - current set of clusters of step n (before combining two clusters)
tree.csv - each line in this file corresponds to two clusters that are going to be combined and MI between them 
"""
import os
import math
import prog4
import common
WORK_DIR = common.WORK_DIR + "exp/"

profile_freq = prog4.get_profile_freq()
profile_pair_freq = prog4.get_profile_pair_freq()
clusters = prog4.get_clusters()
initial_clusters = prog4.get_clusters()


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

def remove_orphans():
    covered_clusters = []
    for i_c in initial_clusters:
        best_MI = -100
        best_final_cluster = set()
        for c in clusters:
            MI = prog4.calc_MI_fast(i_c,c,profile_pair_freq,profile_freq) # calculate MI between initial profile combination (i_c) and each of the current combinations (c)
            if MI > best_MI:
                best_MI = MI
                best_final_cluster = c
        best_final_cluster = sorted(list(best_final_cluster))
        covered_clusters.append(best_final_cluster)

    #now we have to remove duplicates and transform a list of lists into a list of sets
    s = []
    clusters[:] = [] # empty the list
    for i in covered_clusters:
       if i not in s:
          s.append(i)
          clusters.append(set(i))

tree_file = open(WORK_DIR + "tree.csv","w")
step = 1
C = {}  # Mutual information scores for pairs of clusters
while (len(clusters) > 1):
    print "step " + str(step)
    remove_orphans()
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
