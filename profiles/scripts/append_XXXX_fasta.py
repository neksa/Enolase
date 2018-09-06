#!/usr/bin/env python

import string
import sys
import re

input_handle = sys.stdin
output_handle = sys.stdout

sequence = ""
for line in input_handle:
    if line[0:1] == ">":
        if sequence != "":
            output_handle.write("XXXXXXXXXXXX" + sequence.rstrip() + "XXXXXXXXXXXX\n")
        output_handle.write(line)
        sequence = ""
    else:
        sequence += line

if sequence != "":
    output_handle.write("XXXXXXXXXXXX" + sequence.rstrip() + "XXXXXXXXXXXX\n")
