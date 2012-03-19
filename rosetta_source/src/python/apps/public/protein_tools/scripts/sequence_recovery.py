#!/usr/bin/env python2.5

import Bio.PDB
from optparse import OptionParser

from rosettautil.protein import util
from rosettautil.rosetta import rosettaScore

class mutation:
    def __init__(self,resno,nat,des,delta,nat_score,des_score):
        self.native=nat
        self.designed=des
        self.delta=delta
        self.native_score=nat_score
        self.designed_score = des_score
        self.resno = resno

usage = "%prog [options] native.pdb designed.pdb"
parser = OptionParser(usage)
parser.add_option("--best",dest="best",help="print out the top n mutations in terms of total energy, default = 10",default=10)
parser.add_option("--worst",dest="worst",help="print out the bottom n mutations in terms of total energy, default = 10",default=10)
(options,args) = parser.parse_args()

native_struct = util.load_pdb(args[0])
designed_struct = util.load_pdb(args[1])

native_residues = native_struct.get_residues()
designed_residues = designed_struct.get_residues()

native_scores = rosettaScore.ScoreTable(args[0])
designed_scores = rosettaScore.ScoreTable(args[1])



native_sequence = ""
designed_sequence = ""
concensus = ""

conserved_count = 0.0
total_count = 0.0

mutations = []

for native_res, designed_res in zip(native_residues,designed_residues):
    try:
        native_name = Bio.PDB.Polypeptide.three_to_one(native_res.get_resname())
        native_sequence += (native_name)
    except KeyError:
        continue

    try:
        designed_name = Bio.PDB.Polypeptide.three_to_one(designed_res.get_resname())
        designed_sequence += (designed_name)
    except KeyError:
        continue
        
    native_score = native_scores.get_score(native_res.get_id()[1],"total")
    designed_score = designed_scores.get_score(designed_res.get_id()[1],"total")
    delta_score = designed_score-native_score

    if native_name == designed_name:
        concensus += "*"
        conserved_count += 1.0
    else:
        concensus += " "
        mutations.append(mutation(native_res.get_id()[1],native_name,designed_name,delta_score,native_score,designed_score))
    total_count += 1.0

print native_sequence
print designed_sequence
print concensus
print "\n",
print "percent sequence recovery:",(conserved_count/total_count)*100

sorted_mutations = sorted(mutations,key=lambda mut: mut.delta)
best_mutations = sorted_mutations[0:int(options.best)]
worst_mutations = sorted_mutations[len(sorted_mutations)-int(options.worst):len(sorted_mutations)]

print "printing the best",options.best,"mutations"
for mut_record in best_mutations:
    print mut_record.resno,mut_record.native_score, mut_record.native,"->",mut_record.designed,mut_record.designed_score
print "\n",

print "printing the worst",options.worst,"mutations"
for mut_record in worst_mutations:
    print mut_record.resno, mut_record.native_score, mut_record.native,"->",mut_record.designed,mut_record.designed_score
print "\n"
