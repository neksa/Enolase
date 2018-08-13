#!/usr/bin/env python
"""
counts active clusters at each step

"""
import prog4
import prog8
import sys
import common
WORK_DIR = common.WORK_DIR

profile_freq = prog4.get_profile_freq()
profile_pair_freq = prog4.get_profile_pair_freq()
new_file = open(WORK_DIR + "active_clusters.csv","w")
new_file.write("number_of_clusters\tnumber_of_active_clusters\n")
initial_clusters = prog8.find_initial_clusters_and_their_info_from_table2()

for step in range(260,0,-1):
    try:
        current_clusters = prog8.find_clusters_on_a_given_step(step)
        profile_freq = prog4.get_profile_freq()
        profile_pair_freq = prog4.get_profile_pair_freq()

        covered_clusters = set()
        for i_c in initial_clusters:
            best_MI = -100
            best_final_cluster = ""
            for c in current_clusters:
                MI = prog4.calc_MI_fast(set(i_c.split(",")),set(c.split(",")),profile_pair_freq,profile_freq) # calculate MI between initial profile combination (i_c) and each of the current combinations (c)
                if MI > best_MI:
                    best_MI = MI
                    best_final_cluster = c
            covered_clusters.add(best_final_cluster)
        print step, " step"
        new_file.write(str(len(current_clusters)) + '\t' + str(len(covered_clusters)) + "\n")
            
    except IOError:
        print step,'step doesn\'t exist'
    