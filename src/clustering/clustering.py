#!/usr/bin/env python

import math
import os
from common import *


################# DISTANCE BASED ON NORMALIZED POITWISE MUTUAL INFORMATION ######################
def calc_MI(c1, c2, F, G):
    """
    profile combinations - c1 and c2
    F - frequencies of profiles
    G - frequencies of profile pairs
    """
    MI = 0.0
    for p in c1:
        for q in c2:
            if G[p][q] > 0:
                NF = - math.log(G[p][q], 2)
                MI += math.log(G[p][q] / (F[p] * F[q]), 2) / NF

    MI = MI / (len(c1) * len(c2))
    return MI


def calc_MI_matrix(combinations, F, G):
    matrix = {}
    for i, c1 in enumerate(combinations):
        for j, c2 in enumerate(combinations):
            if j >= i:
                continue
            matrix[(i, j)] = calc_MI(c1, c2, F, G)  # round 4
    return matrix


def write_MI_matrix(combinations, C, step):
    if not os.path.exists(WORK_DIR + str(step)):
        os.makedirs(WORK_DIR + str(step))
    with open(WORK_DIR + str(step) + "/MI_matrix.csv", "w") as o:
        for pair, MI in C.iteritems():
            i, j = pair
            c1 = ",".join(map(str, sorted(list(combinations[i]))))
            c2 = ",".join(map(str, sorted(list(combinations[j]))))
            o.write("{}:{}:{}\n".format(c1, c2, MI))


######### AGGLOMERATIVE CLUSTERING ###################
def agglomerative_clustering(F, G, C, no_orphans=False):
    """
    F - profile frequencies
    G - profile pair frequencies
    C - combinations of profiles
    """
    proteins = C[:]
    # singletons = [i for i, c in enumerate(C) if len(c) == 1]

    # remove singltone "C" of profiles
    C = [c for i, c in enumerate(C) if len(c) > 1]

    with open(WORK_DIR + "tree.csv", "w") as o:
        step = 1
        M = {}  # Mutual information scores for pairs of combinations
        while len(C) > 1:
            try:
                for i, c1 in enumerate(C):
                    for j, c2 in enumerate(C):
                        if i > j and c1.issubset(c2) and c2.issubset(c1):
                            del C[j]
                            raise Exception("Merged a pair of identical C: {} {}".format(c1, c2))
            except Exception as e:
                print e
                continue

            if no_orphans:
                C = remove_orphans(proteins, C, F, G)

            print "step ", step
            M = calc_MI_matrix(C, F, G)
            pair, max_MI = max(M.iteritems(), key=lambda x: x[1])
            i, j = pair
            c1 = ",".join(map(str, sorted(list(C[i]))))
            c2 = ",".join(map(str, sorted(list(C[j]))))

            print " argmax", round(max_MI, 3), c1, c2
            o.write("{}:{}:{}\n".format(c1, c2, max_MI))

            write_MI_matrix(C, M, step)
            write_combinations(C, step)
            # merge most similar C
            C.append(C[i] | C[j])
            # remove original C
            C = [cl for k, cl in enumerate(C) if k not in (i, j)]
            step += 1
    # C = calc_MI_matrix(combinations, F, G)
    # create_MI_file_for_this_step(combinations, C, step)
    # create_file_with_information_about_each_cluster_for_this_step(combinaitons, step)


def remove_orphans(proteins, combinations, F, G):
    not_orphans = set()
    for c_prot in proteins:
        best_MI = -10000.0
        best_i = 0
        best_c = set()
        for i, c in enumerate(combinations):
            MI = calc_MI(c_prot, c, F, G)
            if MI > best_MI:
                best_MI = MI
                best_i = i
                best_c = c
            # print MI, c_prot, c
        # print "BEST", c_prot, best_c, best_MI
        not_orphans.add(best_i)
    print "Removed N_orphans", len(combinations) - len(not_orphans)
    return [c for i, c in enumerate(combinations) if i in not_orphans]


######################### CLUSTERING by STEP #####################
def find_steps(l):
    res = []
    list_of_visited_sizes_of_clusters = []

    for number_of_clusters in l:
        number_of_clusters = int(number_of_clusters)
        if number_of_clusters in list_of_visited_sizes_of_clusters:
            continue
        for i in range(1, 329):
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
                print "No step" + str(number_of_clusters)
    res = set(res)
    res = sorted(list(res), reverse=True)
    res = map(str, res)
    return res


def find_initial_clusters_and_their_info_from_table2():
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

    res = {}
    fname = WORK_DIR + "table2.csv"
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Profiles") or line.startswith("Total"):
                continue
            profiles, proteins, proteins_by_subgroups, rest = line.split(";", 3)
            res[profiles] = {}
            res[profiles]['size'] = int(proteins)
            res[profiles]['subgroups'] = translate_subgroups_to_dictionary(proteins_by_subgroups)
            res[profiles]['prediction'] = translate_prediction_to_int(rest)
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
            cl1, cl2, MI = line.split(":")
            MI_sum += float(MI)
            number_of_lines += 1
        average_MI = MI_sum/float(number_of_lines)
    return average_MI


def write_combinations(combinations, step):
    if not os.path.exists(WORK_DIR + str(step)):
        os.makedirs(WORK_DIR + str(step))
    with open(WORK_DIR + str(step) + "/combinations.csv", "w") as o:
        for c in combinations:
            o.write(",".join(map(str, sorted(list(c)))) + '\n')

# def remove_orphans(cc):
#     """
#     cc - current clusters 
#     """
#     initial_clusters = find_initial_clusters_and_their_info_from_table2()
#     current_clusters = cc  #list os sets
#     profile_freq = get_profile_freq()
#     profile_pair_freq = get_profile_pair_freq()
#     covered_clusters = []
#     for i_c in initial_clusters:
#         best_MI = -100
#         best_final_cluster = set()
#         for c in current_clusters:
#             MI = calc_MI(set(i_c.split(",")), c, profile_pair_freq, profile_freq)  # calculate MI between initial profile combination (i_c) and each of the current combinations (c)
#             if MI > best_MI:
#                 best_MI = MI
#                 best_final_cluster = c
#         best_final_cluster = sorted(list(best_final_cluster))
#         covered_clusters.append(best_final_cluster)

#     #now we have to remove duplicates and transform a list of lists into a list of sets
#     s = []
#     new_list = []
#     for i in covered_clusters:
#         if i not in s:
#             s.append(i)
#         new_list.append(set(i))
#     return new_list
###########################################


def count_active_clusters():
    profile_freq = get_profile_freq()
    profile_pair_freq = get_profile_pair_freq()
    with open(WORK_DIR + "active_clusters.csv", "w") as new_file:
        new_file.write("number_of_clusters\tnumber_of_active_clusters\n")
        initial_clusters = prog8.find_initial_clusters_and_their_info_from_table2()

        for step in range(260, 0, -1):
            try:
                current_clusters = find_clusters_on_a_given_step(step)
                profile_freq = get_profile_freq()
                profile_pair_freq = get_profile_pair_freq()

                covered_clusters = set()
                for i_c in initial_clusters:
                    best_MI = -100
                    best_final_cluster = ""
                    for c in current_clusters:
                        # calculate MI between initial profile combination (i_c) and each of the current combinations (c)
                        MI = calc_MI(
                            set(i_c.split(",")),
                            set(c.split(",")),
                            profile_pair_freq,
                            profile_freq)
                        if MI > best_MI:
                            best_MI = MI
                            best_final_cluster = c
                    covered_clusters.add(best_final_cluster)
                print step, " step"
                new_file.write(str(len(current_clusters)) + '\t' + str(len(covered_clusters)) + "\n")

            except IOError:
                print step, 'step doesn\'t exist'


