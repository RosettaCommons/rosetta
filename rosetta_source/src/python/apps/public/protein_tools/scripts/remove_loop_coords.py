#!/usr/bin/env python2.5

import array
import Bio.PDB
from optparse import OptionParser

from rosettautil.protein import util
from rosettautil.rosetta import loops
from rosettautil.util import fileutil

usage = "%prog [options] loopfile.txt input.pdb output.pdb"
parser=OptionParser(usage)
(options,args)=parser.parse_args()

loop_manager = loops.RosettaLoopManager()
loop_manager.read(args[0])
input_struct = util.load_pdb(args[1])

zero_triplet = array.array('f',[0.0,0.0,0.0])

for atom in input_struct.get_atoms():
    resnum = atom.get_parent().get_id()[1]
    if loop_manager.is_res_in_loop(resnum):
        atom.set_coord(zero_triplet)
        atom.set_occupancy(-1.0)

pdb_io = Bio.PDB.PDBIO()
pdb_io.set_structure(input_struct)
outfile = fileutil.universal_open(args[2],'w')
pdb_io.save(outfile)
outfile.close()