#!/usr/bin/env python
#
# filter matrix by id
#

import string
import random
import sys
import math
import re

### parameters ###
if len(sys.argv)-1 < 3:
    print "Usage:       ./profile_analysis.py <input matrix filename> <output matrix filename> <starting_number>"
    sys.exit(1)

matrix_filename = sys.argv[1]
output_filename = sys.argv[2]
starting_id =  int(sys.argv[3])

print "Matrix filename = ", matrix_filename

ID = 0
matrix_counter = starting_id

fin = open(matrix_filename);
fout = open(output_filename, "w");

matrix_buf = ""

for line in fin:
    if line[0:6] == "MATRIX":
        ID = 0
        id_found = False
        pairs = line.split(" ")
        for pair in pairs:
            if pair.find("=") != -1:
                (var, value) = pair.split("=", 1)
                if var == "ID":
                    ID = int(value)
                    id_found = True
        if id_found:
            print "found", ID, "replaced by", matrix_counter
            line = line.replace("ID="+str(ID),"ID="+str(matrix_counter))
        else:
            line = line.replace("MATRIX", "MATRIX ID="+str(matrix_counter))
        ID = matrix_counter

    matrix_buf += line

    if line[0:3] == "END":
        fout.write(matrix_buf)
        matrix_counter += 1
        matrix_buf = ""

fin.close()
fout.close()

print "DONE"
