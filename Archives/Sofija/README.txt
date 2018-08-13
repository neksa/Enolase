table2.csv - is a file with all needed information about occuring initial profile combinations. This is an input file for the clustering.
___________________________________________________________________________
Programs that do not do anything, but include some procedures that are used by other programs

prog8.py
common.py - includes directory path
___________________________________________________________________________
Main programs

1) prog4.py
Should be run first.
Calculates frequecies for each profile as well as profile pairs and generates files: profile_freq.csv and profile_pair_freq.csv
It also includes calc_MI(..) method, which defines how we calculate mutual information between two combinations.
(That is the reason why this program is imported into many other programs.)
calc_MI_fast(..) is similar to calc_MI(..). The difference is that calc_MI_fast(..) runs faster, because it don't need to read a file with frequencies (because frequency-lists are send with parameters).

2) prog1.py

Performs clustering. For each step it creates a folder with step number as a name.
In each step-folder two files are being generated: MI_matrix.csv, clusters.csv(just a list of profile combinations a this step before clustering) 
Also creates a file tree.csv. In this file each line shows which two profile combinations are combined.
___________________________________________________________________________
the rest can be run in a random order

3) prog2.py

Creates file for edges (for cytoscape)

4) prog15.py

Counts active clusters
Output file: active_clusters.csv

5) prog14.py

Calculates FDR
Output file: FDR.csv
