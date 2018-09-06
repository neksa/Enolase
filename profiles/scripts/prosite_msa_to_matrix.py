#!/usr/bin/env python
# PROSITE msa files convert to profile matrix (format 3)
# Alexander Goncearenco 2011
#
# scan directory for *.msa files
# load MSA in fasta format
# set up a matrix in memory
# write sequence logo using Weblogo3 library
# target is 30-residue long profile matrices
# center short profiles in the middle of 30-residue matrix
# remove excessive gaps (>50%) from the matrices
# write variable length PFM frequency profiles format 4
# inform about too long profiles that have to be split

import os
import sys
import glob

import weblogolib
from weblogolib import *
#from weblogolib import  parse_prior, GhostscriptAPI
#from weblogolib.color import *
#from weblogolib.colorscheme import *
#from StringIO import StringIO

from numpy import *
#array, asarray, float64, ones, zeros, int32,all,any, shape

#from corebio import seq_io
from corebio.seq import *
from corebio.matrix import AlphabeticArray

#from subprocess import *
#from pkg_resources import resource_stream

#from corebio.moremath import entropy
#from math import log, sqrt

#print protein_alphabet

profile_alphabet = "ACDEFGHIKLMNPQRSTVWY"

from corebio.matrix import Motif
from corebio.moremath import *

def change_alphabet(m):
    L = len(m)
    new_profile = zeros((L, len(profile_alphabet)), float64)
    for j, amino_acid in enumerate(profile_alphabet):
        j_orig = protein_alphabet.ords(amino_acid)[0]
        new_profile[...,j] = m[...,j_orig]
    return new_profile


def counts_to_PFM(c):
    # normalize on each position to adjust for alphabet conversion losses
    L = len(c)
    PFM = c
    for i in range(L):
        PFM[i] /= sum(PFM[i])
    return PFM

def print_logo(freqs, fname, title, bg_prior=False):
    if bg_prior:
        data = LogoData.from_counts(profile_alphabet, freqs, prior=B)
    else:
        data = LogoData.from_counts(profile_alphabet, freqs, prior=None)
        #data = LogoData.from_counts(profile_alphabet, freqs, prior=ones((len(profile_alphabet)), float64)/len(profile_alphabet))
    #print data

    options = LogoOptions()
    #options.logo_title = "Profile " + str(ID) + " (K=" + str(K) + ")"
    options.logo_title = title
    options.show_yaxis = True
    options.show_errorbars = False
    options.show_title = True
    #options.scale_width = False

    options.stacks_per_line = 50
    options.stack_width = std_sizes['large']
    options.show_fineprint = False
    options.logo_margin = 0.5
    options.creator_text =''

    options.color_scheme = std_color_schemes['chemistry']
    format = LogoFormat(data, options)

    fout = open(fname, 'w')
    png_formatter(data, format, fout)
    fout.close()

"""
PROFILE 3
BEGIN
MATRIX K=45 L=2
   2    A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S    T    V    W    Y    .
   0    1    0    0    2    0    1    1    1    1    0    1    2    0    0    0    0    0    1    0    1   33
   1    1    0    0    0    2    0    0    0    1    1    0    1    2    0    1    0    2    1    0    0   33

PROFILE 4
BEGIN
MATRIX ID=0 K=30 L=30
30        A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
 0 0.089482 0.017268 0.003140 0.009419 0.014129 0.006279 0.006279 0.142857 0.020408 0.031397 0.017268 0.001570 0.004710 0.007849 0.006279 0.084772 0.260597 0.252747 0.001570 0.021978
 1 0.058085 0.014129 0.023548 0.026688 0.001570 0.014129 0.020408 0.144427 0.056515 0.059655 0.009419 0.043956 0.036107 0.010989 0.056515 0.081633 0.054945 0.233909 0.010989 0.042386
"""
def write_matrix(m, K, origin, fname):
    L = len(m)

    f = open(fname, "w")
    f.write("PROFILE 4\n")
    f.write("BEGIN\n")
    f.write("ORIGIN %s\n" % (origin))
    f.write("MATRIX ID=0 K=%d L=%d\n" % (K, L))

    f.write("%2d" % (L))
    for aa in profile_alphabet:
        f.write(" %8s" % (aa))
    f.write("\n")

    for i in range(L):
        f.write("%2d" % (i))
        for j, aa in enumerate(profile_alphabet):
            f.write(" %1.6f" % (m[i,j]))
        f.write("\n")

    f.write("END\n")
    f.close()

def read_bg_composition(fname):

    B = zeros(len(profile_alphabet), float64)

    if os.path.exists(fname):
        fc = open(fname)
        for line in fc:
            (aa, bg) = line.split(",")
            bg = float(bg)/100.0
            B[profile_alphabet.index(aa)] = float(bg)
        fc.close()
    else:
        # Attention: Archaeal frequencies by default
        B[profile_alphabet.index('A')] = 0.07846065
        B[profile_alphabet.index('R')] = 0.05395213
        B[profile_alphabet.index('N')] = 0.03854661
        B[profile_alphabet.index('D')] = 0.05671408
        B[profile_alphabet.index('C')] = 0.00983962
        B[profile_alphabet.index('E')] = 0.07457194
        B[profile_alphabet.index('Q')] = 0.02283004
        B[profile_alphabet.index('G')] = 0.07429558
        B[profile_alphabet.index('H')] = 0.01708365
        B[profile_alphabet.index('I')] = 0.07471997
        B[profile_alphabet.index('L')] = 0.09529720
        B[profile_alphabet.index('K')] = 0.05845627
        B[profile_alphabet.index('M')] = 0.02372575
        B[profile_alphabet.index('F')] = 0.03902878
        B[profile_alphabet.index('P')] = 0.04283092
        B[profile_alphabet.index('S')] = 0.06101052
        B[profile_alphabet.index('T')] = 0.05260790
        B[profile_alphabet.index('W')] = 0.01027624
        B[profile_alphabet.index('Y')] = 0.03727149
        B[profile_alphabet.index('V')] = 0.07847484
        # normalize
        B /= sum(B)
    return B

"""
GR      AC      DO      ID      LN      TY      PA
ATP     PS00017 PDOC00017       ATP_GTP_A       0       PATTERN [AG]-x(4)-G-K-[ST]
"""
def read_prosite_map_groups(fname):
    f = open(fname)
    first_line = True
    prosite_map = {}
    for line in f:
        if first_line:
            first_line = False
            continue
        (GR, AC, DO, ID, LN, TY, PA) = line.split("\t")
        prosite_map[AC] = {'GR':GR, 'AC':AC, 'DO':DO, 'ID':ID, 'LN':LN, 'TY':TY, 'PA':PA.strip()}
        #print AC, prosite_map[AC]
    f.close()
    return prosite_map


def adjust_matrix(PFM, B, L, M):
    adjusted_matrix = zeros((M, len(profile_alphabet)), float64)
    spacer = (M-L)/2
    # left background spacer
    for i in range(spacer):
        adjusted_matrix[i] = B
    # copy values from the original matrix in the middle of the new one
    for i in range(L):
        adjusted_matrix[spacer + i] = PFM[i]
    # right background spacer
    for i in range(spacer+i, max(M, L + 2*spacer)):
        adjusted_matrix[i] = B

    return adjusted_matrix

##################################################
#version = "20.75"
#database_path = "/net/workspace/users/agoncear/Prosite/"+version+"/"

B = read_bg_composition("composition.csv")
#print B

print "Reading Prosite PS-PDOC map with groups"
prosite_map = read_prosite_map_groups("PS_map_groups.tab")

print "Reading alignments"
ps_length = {}

for aln_file in glob.glob("selected/alignments/PS*.msa"):
#for aln_file in glob.glob("selected/alignments/PS00112.msa"):
    f = open(aln_file)
    AC = aln_file[aln_file.find('PS'): aln_file.find('.msa')]
    #print "reading", aln_file
    sequence = ""
    freqs = zeros((100, len(protein_alphabet))) # just to init
    K = 0
    L = 0
    for line in f:
        if line[0:1] == ">":
            if sequence != "":
                K += 1
                sequence = sequence.upper().replace(".", "-")
                L = len(sequence)
                for i in range(0, L):
                    #print i, sequence[i].upper()
                    freqs[i, protein_alphabet.ords(sequence[i])[0]] += 1

                sequence = ""
        else:
            sequence += line.strip()
    f.close()
    #print AC, prosite_map[AC]['GR'], prosite_map[AC]['DO'], K, L
    #print K

    #rtrim the matrix
    freqs = delete(freqs, s_[L::], 0)

    positions_to_remove = 0
    trimmed_freqs = freqs
    for i in reversed(range(L)):
        N_gaps = freqs[i, protein_alphabet.ords('-')[0]]
        proportion_gaps = N_gaps / K
        if proportion_gaps > 0.5:
            #print AC, i, proportion_gaps
            positions_to_remove += 1
            trimmed_freqs = delete(trimmed_freqs, i, 0);
        #for a in protein_alphabet:
        #print i, a, freqs[i, protein_alphabet.ords(a)[0]]


    origin = "%s; %s; %s; 1" % (prosite_map[AC]['GR'], prosite_map[AC]['DO'], AC)

    profile = change_alphabet(trimmed_freqs)
    PFM = counts_to_PFM(profile)
    #print_logo(PFM, "selected/3_after_trunc_logo/"+AC+".png", origin)

    original_L = L
    L = original_L - positions_to_remove

    if positions_to_remove > 0:
        print prosite_map[AC]['GR'], prosite_map[AC]['DO'], AC, L, positions_to_remove

    #print PFM
    if L > 50:
        print "Attention: the profile %s is longer than 50 residues (%d)" % (AC, L)
        adjusted_matrix = PFM
    elif L > 30:
        # target L = 50
        adjusted_matrix = adjust_matrix(PFM, B, L, 50)
    else:
        # target L = 30
        adjusted_matrix = adjust_matrix(PFM, B, L, 30)

    #print and save adjusted_matrix
    print_logo(adjusted_matrix, "selected/4_after_cut_logo/"+AC+".png", origin)
    write_matrix(adjusted_matrix, K, origin, "selected/4_converted_matrices/"+AC+".matrix")
