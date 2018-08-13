#!/usr/bin/env python
"""
USEFUL METHODS
translates number of clusters into numbers of steps using clusters.csv-files
and other useful methods
"""
import common
WORK_DIR = common.WORK_DIR

def find_steps(l):
    res = []
    list_of_visited_sizes_of_clusters = []

    for number_of_clusters in l:
        number_of_clusters = int(number_of_clusters)
        if number_of_clusters in list_of_visited_sizes_of_clusters: continue
        for i in range(1,329):
            try:
                lines_number = 0
                with open(WORK_DIR + str(i) + '/clusters.csv') as infp:
                    for line in infp:
                        if line.strip():
                            lines_number += 1
                if lines_number == number_of_clusters:
                    res.append(i)
                    list_of_visited_sizes_of_clusters.append(number_of_clusters)
                    break
                if lines_number < number_of_clusters:
                    res.append(i)
                    list_of_visited_sizes_of_clusters.append(lines_number)
                    break
            except IOError:
                print "prog8.py: No step" + str(number_of_clusters)
    res = set(res)
    res = sorted(list(res), reverse=True)
    res = map(str, res)
    return res

def find_initial_clusters_and_their_info_from_table2():
    res = {}
    fname = WORK_DIR + "table2.csv"
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Profiles"): continue
            if line.startswith("Total"): continue
            profiles,proteins,proteins_by_subgroups,rest = line.split(";",3)
            res[profiles]={}
            res[profiles]['size'] = int(proteins)
            res[profiles]['subgroups'] = translate_subgroups_to_dictionary(proteins_by_subgroups)
            res[profiles]['prediction'] = translate_prediction_to_int(rest)
    return res

def translate_prediction_to_int(l):
    return int(l.split("(")[0]) + int(l.split("(")[1].split(")")[0])

def translate_subgroups_to_dictionary(l):
    res = {}
    if l == "":
        return res
    else:
        l = l.split("_")
        for gr in l:
            subgroup, proteins_in_subgroup = gr.split(":")
            res[subgroup] = int(proteins_in_subgroup)
        return res

def find_clusters_on_a_given_step(s):
    res = []
    fname = WORK_DIR + str(s) + "/clusters.csv"
    with open(fname) as f:
        for line in f:
            line = line.strip()
            res.append(line)
    return res

def find_MI_mean_on_a_given_step(s):
    MI_sum = 0
    number_of_lines = 0
    fname = WORK_DIR + str(s) + "/MI_matrix.csv"
    with open(fname) as f:
        for line in f:
            line = line.strip()
            cl1,cl2,MI = line.split(":")
            MI_sum += float(MI)
            number_of_lines += 1
        average_MI = MI_sum/float(number_of_lines)
    return average_MI

def remove_orphans(cc): ## cc - current clusters 
    initial_clusters = find_initial_clusters_and_their_info_from_table2()
    current_clusters = cc #list os sets
    import prog4
    profile_freq = prog4.get_profile_freq()
    profile_pair_freq = prog4.get_profile_pair_freq()
    covered_clusters = []
    for i_c in initial_clusters:
        best_MI = -100
        best_final_cluster = set()
        for c in current_clusters:
            MI = prog4.calc_MI_fast(set(i_c.split(",")),c,profile_pair_freq,profile_freq) # calculate MI between initial profile combination (i_c) and each of the current combinations (c)
            if MI > best_MI:
                best_MI = MI
                best_final_cluster = c
        best_final_cluster = sorted(list(best_final_cluster))
        covered_clusters.append(best_final_cluster)

    #now we have to remove duplicates and transform a list of lists into a list of sets
    s = []
    new_list = []
    for i in covered_clusters:
       if i not in s:
          s.append(i)
          new_list.append(set(i))
    return new_list