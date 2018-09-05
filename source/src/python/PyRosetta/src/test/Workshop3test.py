#! /usr/bin/python
# List of commands used in PyRosetts Workshop #3

from __future__ import print_function

from math import *

from rosetta import *
from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.teaching import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

# Scoring Poses
ras = pose_from_file("../test/data/workshops/6Q21.clean.pdb")

scorefxn = create_score_function("ref2015")

print( scorefxn )

scorefxn2 = ScoreFunction()
scorefxn2.set_weight(fa_atr, 1.0)
scorefxn2.set_weight(fa_rep, 1.0)

print( scorefxn(ras) )

scorefxn.show(ras)

print( ras.energies().show(24) )

r1 = ras.residue(24)
r2 = ras.residue(20)
a1 = r1.atom_index("N")
a2 = r2.atom_index("O")


etable_atom_pair_energies(r1, a1, r2, a2, scorefxn)


hbond_set = hbonds.HBondSet()
ras.update_residue_neighbors()
hbonds.fill_hbond_set(ras, False, hbond_set)
hbond_set.show(ras)

hset = ras.get_hbonds()
hset.show(ras)
hset.show(ras, 24)

pose = pose_from_file("../test/data/workshops/1YY9.clean.pdb")

rsd1_num = pose.pdb_info().pdb2pose('D', 102)
rsd2_num = pose.pdb_info().pdb2pose('A', 408)
print( rsd1_num )
print( rsd2_num )
rsd1 = pose.residue(rsd1_num)
rsd2 = pose.residue(rsd2_num)

emap = EMapVector()
scorefxn.eval_ci_2b(rsd1, rsd2, pose, emap)
print( emap[fa_atr] )
print( emap[fa_rep] )
print( emap[fa_sol] )

pymol = PyMOLMover()
ras.pdb_info().name("ras")
pymol.send_energy(ras)
pymol.send_energy(ras, "fa_atr")
pymol.send_energy(ras, "fa_sol")
pymol.update_energy(True)
pymol.energy_type(fa_atr)
pymol.apply(ras)
# no longer supported: pymol.label_energy(ras, "fa_rep")
# no longer supported: pymol.send_hbonds(ras)
