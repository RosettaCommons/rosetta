#! /usr/bin/python
# List of commands used in PyRosetts Workshop #3

from math import *

from rosetta import *
init()

# Scoring Poses
ras = pose_from_pdb("test/data/workshops/6Q21.clean.pdb")

scorefxn = create_score_function("standard")

print scorefxn

scorefxn2 = ScoreFunction()
scorefxn2.set_weight(fa_atr, 1.0)
scorefxn2.set_weight(fa_rep, 1.0)

print scorefxn(ras)

scorefxn.show(ras)

print ras.energies().show(24)

hset = get_hbonds(ras)
hset.show(ras)
hset.show(ras, 24)

pose = pose_from_pdb("test/data/workshops/1YY9.clean.pdb")

rsd1_num = pose.pdb_info().pdb2pose('D', 102)
rsd2_num = pose.pdb_info().pdb2pose('A', 408)

rsd1 = pose.residue(rsd1_num)
rsd2 = pose.residue(rsd2_num)

emap = EMapVector()
scorefxn.eval_ci_2b(rsd1, rsd2, pose, emap)
print emap[fa_atr]
print emap[fa_rep]
print emap[fa_sol]

