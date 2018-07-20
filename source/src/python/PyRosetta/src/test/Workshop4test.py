#! /usr/bin/python
# :noTabs=true:
# List of commands used in PyRosetts Workshop #4

from __future__ import print_function

import sys
#if sys.platform == "darwin": sys.exit(0)  # skipping this test on Mac OS due to memory error (*** error: can't allocate region)


# A simple de Novo Folding Algorithm
# import modules
import random, math

from rosetta import *
from pyrosetta import *
from pyrosetta.teaching import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

# constants
PHI = 0
PSI = 1
kT = 1.0

# subroutines/functions

def random_move(pose):
    # select random residue
    res = random.randint(1, pose.total_residue())

    # select and set random torsion angle distributed around old angle
    if random.randint(PHI, PSI) == PHI:
        torsion = pose.phi(res)
        a = random.gauss(torsion, 25)
        pose.set_phi(res, a)
    else:
        torsion = pose.psi(res)
        a = random.gauss(torsion, 25)
        pose.set_psi(res, a)

# initialize pose objects
last_pose = Pose()
low_pose = Pose()

# create poly-A chain and set all peptide bonds to trans
p = pose_from_sequence("AAAAAAAAAA", "fa_standard")
for res in range(1, p.total_residue() + 1):
    p.set_omega(res, 180)

# use the PyMOLMover to echo this structure to PyMOL
pmm = PyMOLMover()
pmm.apply(p)

# set score function to include Van der Wals and H-bonds only
score = ScoreFunction()
score.set_weight(fa_atr, 0.8)
score.set_weight(fa_rep, 0.44)
score.set_weight(hbond_sr_bb, 1.17)

# initialize low score objects
low_pose.assign(p)
low_score = score(p)

for i in range(100):
    last_score = score(p)
    last_pose.assign(p)

    random_move(p)

    new_score = score(p)
    print( "Iteration:", i, "Score:", new_score, "Low Score:", low_score )
    deltaE = new_score - last_score

    # accept if new energy score is improved
    # reject or accept if new energy score is worse based on Metropolis criteria
    if deltaE > 0:
        P = math.exp(-deltaE/kT)  # probability of accepting move diminishes
                             # exponetially with increasing energy
        roll = random.uniform(0.0, 1.0)
        if roll >= P:
            p.assign(last_pose)  # reject pose and reassign previous
            continue

    # if new pose is accepted, store lowest score and associated pose
    if new_score < low_score:
        low_score = new_score
        low_pose.assign(p)

# output files
p.dump_pdb("poly-A_final.pdb")
low_pose.dump_pdb("poly-A_low.pdb")


# Low-Resolution (Centroid) Scoring
ras = pose_from_file("../test/data/workshops/6Q21.clean.pdb")
score2 = get_score_function()
print( score2(ras) )
print( ras.residue(5) )

switch = SwitchResidueTypeSetMover("centroid")
switch.apply(ras)
print( ras.residue(5) )

score3 = create_score_function("score3")
print( score3(ras) )

switch2 = SwitchResidueTypeSetMover("fa_standard")
switch2.apply(ras)
print( ras.residue(5) )

# Protein Fragments
fragset = core.fragment.ConstantLengthFragSet(3)
fragset.read_fragment_file("../test/data/workshops/aat000_03_05.200_v1_3")

movemap = MoveMap()
movemap.set_bb(True)
mover_3mer = protocols.simple_moves.ClassicFragmentMover(fragset, movemap)

pose = pose_from_sequence("RFPMMSTFKVLLCGAVLSRIDAG", "centroid")
for res in range(1, p.total_residue() + 1):
    pose.set_omega(res, 180)

mover_3mer.apply(pose)
