#!/usr/bin/env python2.5

from optparse import OptionParser
from rosettautil.rosetta import rosettaScore
from rosettautil.util import fileutil

usage = "%prog [options] --term=scoreterm --percent=10 --mode=silent input_file"
parser=OptionParser(usage)
parser.add_option("--term",dest="term",help="scoreterm to use",default="")
parser.add_option("--percent",dest="percent",help="percent of structures to output",default=10)
parser.add_option("--mode",dest="mode",help="file input mode, use 'pdb' for list of pdb file paths or 'silent' for a silent file",default="silent")
(options,args) = parser.parse_args()

if options.term=="":
    parser.error("you must specify a scoring term with the --term flag")
if len(args) != 1:
    parser.error("you must specify a silent file or a file containing a list of pdb file paths")

if options.mode=="silent" or options.mode=="Silent":
    scores = rosettaScore.SilentScoreTable()
    scores.add_file(args[0])

    structure_count = len(scores)
    percent = float(options.percent)/100.0
    structs_to_print =int(percent*structure_count)

    sorted_scores = scores.sorted_score_generator(options.term)

    count = 0
    while count < structs_to_print:
        (tag,score) = sorted_scores.next()
        print tag,score
        count += 1
elif options.mode =="pdb" or options.mode=="PDB":
    file_scores = []
    pathfile = fileutil.universal_open(args[0],"r")
    for path in pathfile:
        path = path.rstrip()
        scores = rosettaScore.ScoreTable(path)
        total_score = scores.get_score(0,options.term)
        file_scores.append( (path,total_score) )
    file_scores = sorted(file_scores,key=lambda x: x[1])
    structure_count = len(file_scores)
    percent = float(options.percent)/100.0
    structs_to_print = int(percent*structure_count)
    count = 0
    for i in range(structs_to_print):
        (tag,score) = file_scores[i]
        print tag,score
        
else:
    parser.error("you must specify either 'pdb' or 'silent' with the --mode flag")