#!/usr/bin/env python
"""
creates matrices for profile_freq, profile_pair_freq and list of initial clusters
+ def calc_MI
"""
import common
WORK_DIR = common.WORK_DIR
import math     

def calculate_frequencies():
    fname = WORK_DIR + "table2.csv"

    def weight(N, cluster):
        if N > len(cluster):
            return 1
        return 0

    profile_freq = {}
    pair_freq = {}
    norm_factor = 0.0
    clusters = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Prof"): continue
            if line.startswith("Total"): continue
            cl, proteins, the_rest = line.split(";", 2)
            proteins = int(proteins)
            cluster = set()
            for p in cl.split(","):
                p.replace(" ","")
                cluster |= set([p])

            norm_factor += weight(proteins, cluster)
            for p in cluster:
                if p not in profile_freq:
                    profile_freq[p] = 0.0
                profile_freq[p] += weight(proteins, cluster)
            
            for p in cluster:
                for q in cluster:
                    pair = (p, q)
                    if pair not in pair_freq:
                        pair_freq[pair] = 0.0
                    pair_freq[pair] += weight(proteins, cluster)
            clusters.append(cluster)


    #create file with profile frequencies
    new_file = open(WORK_DIR + "profile_freq.csv","w")
    for p in profile_freq:
        profile_freq[p] /= norm_factor
        new_file.write("%s:%.15f\n" % (p,profile_freq[p]))
    new_file.close()

    profile_pair_freq = {}
    #create file with profiles pair frequencies
    new_file = open(WORK_DIR + "profile_pair_freq.csv","w")
    for p in profile_freq:
        profile_pair_freq[p] = {}
        for q in profile_freq:
            profile_pair_freq[p][q] = 0.0
            if (p,q) in pair_freq:
                profile_pair_freq[p][q] = pair_freq[(p,q)] / norm_factor
            new_file.write("%s:%s:%.15f\n" % (p,q,profile_pair_freq[p][q]))
    new_file.close()
    

if __name__ == '__main__':
    calculate_frequencies()


def get_profile_freq():
    profile_freq = {}
    new_file = open(WORK_DIR + "profile_freq.csv","r")
    for line in new_file:
        p,value = line.split(":")
        profile_freq[p] = float(value)
    new_file.close()
    return profile_freq

def get_profile_pair_freq():
    profile_pair_freq = {}
    new_file = open(WORK_DIR + "profile_pair_freq.csv","r")
    for line in new_file:
        p1,p2,value = line.split(":",2)
        if p1 not in profile_pair_freq:
            profile_pair_freq[p1] = {}
        profile_pair_freq[p1][p2] = float(value)
    new_file.close()
    return profile_pair_freq

def get_clusters():
    clusters = []
    fname = WORK_DIR + "table2.csv"
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Profiles"): continue
            if line.startswith("Total"): continue
            profiles,rest = line.split(";",1)
            cluster = set()
            for p in profiles.split(","):
                p.replace(" ","")
                cluster |= set([p])
            clusters.append(cluster)            
    return clusters

def calc_MI_fast(c1, c2,ppf,pf):
    profile_pair_freq = ppf
    profile_freq = pf
    MI = 0.0
    for p in c1:
        for q in c2:
            if profile_pair_freq[p][q] > 0:

                NF = - math.log(profile_pair_freq[p][q],2)
                MI += math.log( profile_pair_freq[p][q] / (profile_freq[p] * profile_freq[q]), 2) / NF
                
    MI = MI / (len(c1) * len(c2))
    return MI

def calc_MI(c1, c2):
    profile_pair_freq = get_profile_pair_freq()
    profile_freq = get_profile_freq()
    MI = 0.0
    for p in c1:
        for q in c2:
            if profile_pair_freq[p][q] > 0:

                NF = - math.log(profile_pair_freq[p][q],2)
                MI += math.log( profile_pair_freq[p][q] / (profile_freq[p] * profile_freq[q]), 2) / NF
                
    MI = MI / (len(c1) * len(c2))
    return MI
