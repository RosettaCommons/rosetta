#!/usr/bin/env python2.5
from optparse import OptionParser
from rosettautil.protein import util
from rosettautil.protein import pdbStat
from rosettautil.rosetta import rosettaScore
from rosettautil.util import fileutil
import Bio.PDB


def make_table(name_score_rmsd,out_name):
    table = fileutil.universal_open(out_name,'w')
    table.write("file\tscore\tRMSD\n")
    for point in name_score_rmsd:
        table.write(str(point[0])+"\t"+str(point[2])+"\t"+str(point[1])+"\n")
    table.close()

usage = "%prog [options] --native=native.pdb --table=out_table.txt --term=scoreterm list of pdb files"
parser = OptionParser(usage)
parser.add_option("--native",dest="native",help="native structure")
parser.add_option("--CA",dest="ca",help="calculate CA rmsd",default=False,action="store_true")
parser.add_option("--res",dest="residues",help="list of residues to use",default="")
parser.add_option("--table",dest="table",help="output a tab delimited table",default="")
parser.add_option("--chain",dest="chain",help="calculate RMSD of only 1 chain",default="")
parser.add_option("--term",dest="term",help="score term to use")
(options,args)=parser.parse_args()

if len(args) <1:
    parser.error("specify at least 1 protein to compare to native")

print "calculating RMSD between ", len(args), " decoys and the native structure ", options.native
if options.ca:
    print "calculating CA RMSD"
else:
    print "calculating all-atom RMSD"

native = util.load_pdb(options.native)
score_rmsd = []
name_score_rmsd =[]
for decoy_file in args:
    tag = decoy_file.split("/").pop().split(".")[0]
    decoy = util.load_pdb(decoy_file.rstrip())
    score_table = rosettaScore.ScoreTable(decoy_file)
    rms =  pdbStat.calculate_rms(native,decoy,options.ca,options.residues,"",options.chain)
    score = score_table.get_score(0,options.term) #0 is the pose energy
    score_rmsd.append( (rms,score) )
    name_score_rmsd.append( (decoy_file,rms,score) )
    
make_table(name_score_rmsd,options.table)
