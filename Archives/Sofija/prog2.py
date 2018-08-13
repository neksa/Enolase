#!/usr/bin/env python
"""
creates edge file for cytoscape

input:
1/MI_matrix.csv - MI matrix of step 1 (before combining two clusters)

output:
1/cyt_edge_file.csv - edge file for step 1
"""
import common
WORK_DIR = common.WORK_DIR

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

output_file = open(WORK_DIR + "1/cyt_edge_file.csv","w")
output_file.write("cluster1\tIC_IC\tcluster2\tMI-score\n")
with open(fname) as f:
    for line in f:
        line = line.strip()
        cl1, cl2, mi = line.split(":", 2)
        if float(mi) > average_mi:
            output_file.write("%s\tCL_CL\t%s\t%s\n" % (cl1,cl2,mi))
output_file.close()
