#!/usr/bin/env python2.5
import sys
from rosettautil.rosetta import rosettaScore
from rosettautil.protein import util
from rosettautil.util import fileutil
from Bio.PDB import *
from optparse import OptionParser

usage = "%prog input.pdb output.pdb"
parser= OptionParser(usage)
parser.add_option("-n",dest="start",help="residue number to start with, default is 1",default=1)
parser.add_option("--preserve",dest="preserve",help="preserve insertion code and heteroflags",default=False, action="store_true")
parser.add_option("--norestart",dest="norestart",help="don't start renumbering at each chain, default=False",default=False, action="store_true")
parser.add_option("--keep-table",dest="table",help="Preserve the rosetta score table at the bottom of the pdb",action="store_true",default=False)
(options, args) = parser.parse_args()


struct = util.load_pdb(args[0])
try:
    residue_id = int(options.start)
except ValueError:
    sys.exit("residue number specified with -n must be an integer")
    
chain_id = ""
for residue in struct.get_residues():
    chain = residue.get_parent()
    if(chain_id != chain.get_id() and not options.norestart):
        chain_id = chain.get_id()
        residue_id=int(options.start)
    #print chain.get_id()
    if(options.preserve):
        hetero = residue.id[0]
        insert = residue.id[2]
        residue.id=(hetero,residue_id,insert)
    else:
        residue.id=(' ',residue_id,' ')
    residue_id +=1
    
    
io=PDBIO()
io.set_structure(struct)
outfile = fileutil.universal_open(args[1],'w')
io.save(outfile)

if(options.table):
    raw_table = rosettaScore.get_table(args[0])
    #outfile = fileutil.universal_open(args[1],'a')
    outfile.writelines(raw_table)
outfile.close()
