#!/usr/bin/python
##
## make boinc submit files for homologs
##
from sys import argv,stdout,exit,stderr
from os import popen,system
from os.path import dirname, basename, abspath, exists
import string
from glob import glob
from math import sin, cos, degrees, radians, atan
import math

if len(argv) == 1:
    print "EXTRACT_BY_SIDECHAIN_SASA_SCORE"
    print
    print "./extract_by_sidechain_sasa_score.py <sc_sasa_scorefile> <outfile> <cutoff>"
    print
    print "creates a new file, <outfile>.sasa.out, that contains only the structures in outfile which have a sidechain sasa score <= to the cutoff. The default cutoff is zero."
    print
    print "make <sc_sasa_scorefile> using the following command line:"
    print
    print "score_sidechain_sasa.linuxgccrelease  -database <mini_db> -in:file:fullatom -in:file:silent <outfile> -scorefile <sc_sasa_scorefile> -score::sidechain_exposed <list>  -score::sidechain_buried <list>"
    print
    print "( <list> means just list out the desired residue numbers, ie -sidechain_exposed 1 2 3 5 7 9 )"
    print
    exit()

scorefile = open(argv[1]).readlines()
outfile = open(argv[2])
extract = open(argv[2]+".sasa.out", 'w')

cutoff = 0
if len(argv) > 3:
    cutoff = float(argv[3])

tags = {}
for i in scorefile:
    line = i.split()

    if line[1] != "sasa_score":

        sasa_score = float(line[1])
        tag = line[2]
        if sasa_score <= cutoff:
            tags[tag] = 0

l = outfile.readline()
while l != "":
    line = l.split()

    if len(line) > 1:
        if line[0] == "SEQUENCE:" or line[1] == "score":
            extract.write(l)
        else:
            tag = line[len(line)-1]

            if tags.has_key(tag):
                extract.write(l)

    l = outfile.readline()
