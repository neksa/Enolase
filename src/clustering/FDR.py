#!/usr/bin/env python
"""
FDR calculation

1) you don't need to understand what is happening before the second for-loop
2) for each current cluster identify the "true" subgroup, which is the subgroup of the nearest protein(s)
3) for each initial cluster, identify its nearest current cluster
4) count the number of proteins whose subgroup don't match the "true" subgroup of the nearest current cluster
5) divide this number by the total number of the proteins, whose subgroup is identified
"""

from clustering import *
from commmon import *

new_file = open(WORK_DIR + "FDR.csv", "w")
new_file.write("clusters\tFDR\n")

initial_clusters = find_initial_clusters_and_their_info_from_table2()

profile_freq = get_profile_freq()
profile_pair_freq = get_profile_pair_freq()

#this is the second for-loop
for step in range(260, 0, -1):
    try:
        print step, " step"
        current_clusters = prog8.find_clusters_on_a_given_step(step)
        new_file.write(str(len(current_clusters)) + "\t")

        # here I will store the name of the "nearest" subgroup for each of the current clusters.
        # For instance, in figure 6A (https://www.dropbox.com/s/bew8l8vszgkgjs5/fig6.png),
        # p3 is the nearest for C1, and p4 is the nearest for C2.
        best_subgroup_of_final_cluster = {}

        for c in current_clusters:
            # best MI so far for this current cluster
            best_MI = -100
            best_subgroups_in_current_clusters = {}

            for i_c in initial_clusters:
                # calculate MI between an initial cluster and a current cluster
                MI = prog4.calc_MI_fast(set(i_c.split(",")), set(c.split(",")), profile_pair_freq, profile_freq)
                # if MI is higher than the best MI so far, and there exist any proteins with known subgroup in this initial cluster
                if MI > best_MI and initial_clusters[i_c]['subgroups'] != {}:
                    best_MI = MI
                    best_subgroups_in_current_clusters = {}

                # if MI is higher than the best MI so far, and there exist any proteins with known subgroup in this initial cluster
                if MI == best_MI:
                    # keys are the names of the subgroups, values are the numbers of proteins from the given subgroup
                    for subgr in initial_clusters[i_c]['subgroups'].keys():
                        if subgr not in best_subgroups_in_current_clusters:
                            best_subgroups_in_current_clusters[subgr] = 0
                        best_subgroups_in_current_clusters[subgr] += initial_clusters[i_c]['subgroups'][subgr]

            best_subgr = max(best_subgroups_in_current_clusters.keys(), key=lambda x: best_subgroups_in_current_clusters[x])
            # store the name of the subgroup of the proteins which are the nearest to the current cluster.
            # This subgroup is the most represented one among the nearest proteins with known subgroup.
            best_subgroup_of_final_cluster[c] = best_subgr

        # now we know the nearest subgroup of each current cluster,
        # so let's go through initial clusters and find THEIR nearest current cluster.
        # And calculate FP, and then FDR for this particular step.
        total_FP = 0
        total = 0

        for i_c in initial_clusters:
            best_MI = -100
            # nearest current cluster for this initial cluster
            best_final_cluster = ""
            for c in current_clusters:
                MI = prog4.calc_MI_fast(set(i_c.split(",")), set(c.split(",")), profile_pair_freq, profile_freq)
                if MI > best_MI:
                    best_MI = MI
                    best_final_cluster = c

            # now, there are probably proteins from different subgroups in this initial cluster
            for subgr in initial_clusters[i_c]['subgroups'].keys():
                 # but we are only interested in False subgroups (not the nearest subgroup of the nearest current cluster)
                if subgr != best_subgroup_of_final_cluster[best_final_cluster]:
                    total_FP += initial_clusters[i_c]['subgroups'][subgr]
            # here, we count all the proteins, whose subgroup is known
            total += sum(initial_clusters[i_c]['subgroups'].values())

        # now we can calculate FDR

        new_file.write(str(total_FP/float(total)) + "\n")

    except IOError:
        print step, 'step doesnt exist'
