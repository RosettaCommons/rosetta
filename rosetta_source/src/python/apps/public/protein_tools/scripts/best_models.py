#!/usr/bin/env python2.5

from optparse import OptionParser
from rosettautil.rosetta import rosettaScore

usage = "%prog [options] --term=scoreterm silent files"
parser=OptionParser(usage)
parser.add_option("--term",dest="term",help="score term to use")
(options,args) = parser.parse_args()

if len(args) < 1:
    parser.error("you must specify at least one silent file")
#score_gen = scores.score_generator(options.term)
best_models = {} # key is a structure ID, value is a pair in form (tag,score)
for silent_file in args:
    #file = silent_file
    scores=rosettaScore.SilentScoreTable()
    scores.add_file(silent_file)
    score_gen = scores.score_generator(options.term)
    for tag,score in score_gen:
        split_tag = tag.split("_")
        model_id = "_".join(split_tag[0:len(split_tag)-1])
        #file = scores.get_file_from_tag(tag)
        try:
            (current_file,current_best_tag,current_best_score) = best_models[model_id]
        except KeyError:
            best_models[model_id] = (silent_file,tag,score)
            continue
        if score < current_best_score:
            #print "changed"
            best_models[model_id] = (silent_file,tag,score)
            #print best_models
        #print silent_file
        #print file,score , current_best_score
    

print "file","tag",options.term
for tag in best_models:
    print best_models[tag][0],best_models[tag][1],best_models[tag][2]
    
