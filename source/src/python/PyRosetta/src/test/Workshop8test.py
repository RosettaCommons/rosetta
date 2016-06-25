#! /usr/bin/python
# List of commands used in PyRosetts Workshop #2

from __future__ import print_function

from rosetta import *
from pyrosetta import *
from pyrosetta.teaching import *

from rosetta.protocols.loops.loop_closure.ccd import *
from rosetta.protocols.loops.loop_mover.refine import *
from rosetta.protocols.loops.loop_mover.perturb import *
from rosetta.protocols.loops.loop_closure.kinematic_closure import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

pose = pose_from_file("../test/data/test_in.pdb")

# Fold Tree
ft = FoldTree()
ft.add_edge(1, 13, -1)
ft.add_edge(13, 19, -1)
ft.add_edge(13, 26, 1)
ft.add_edge(26, 20, -1)
ft.add_edge(26, 116, -1)

print( ft )
ft.check_fold_tree()

pose.fold_tree(ft)

for res in (10, 13, 16, 23, 26, 30):
    pose.set_phi(res, 180)
    pose.dump_pdb("loop" + str(res) + ".pdb")

pmm = PyMolMover()
pmm.apply(pose)
# no longer supported: pmm.send_foldtree(pose)
# no longer supported: pmm.view_foldtree_diagram(pose, ft)

ft.clear()
ft.simple_tree(116)
ft.new_jump(76, 85, 80)

# Cyclic Coordination Descent (CCD) Loop Closure
movemap = MoveMap()
movemap.set_bb(True)
movemap.set_chi(True)

loop1 = protocols.loops.Loop(15, 24, 19)
add_single_cutpoint_variant(pose, loop1)
ccd = protocols.loops.loop_closure.ccd.CCDLoopClosureMover(loop1, movemap)

ccd.apply(pose)
set_single_loop_fold_tree(pose, loop1)

# Multiple Loops
loop2 = protocols.loops.Loop(78, 83, 80)

loops = protocols.loops.Loops()
loops.add_loop(loop1)
loops.add_loop(loop2)

# Loop Building
reference_pose = pose_from_file("../test/data/test_in.pdb")
score = get_score_function()

import tempfile
output = tempfile.mkstemp()[1]

jd = PyJobDistributor(output, 1, score)

lrms = loop_rmsd(pose, reference_pose, loops, True)
jd.additional_decoy_info = " LRMSD: " + str(lrms)

# High-Resolution Loop Protocol
loop_refine = LoopMover_Refine_CCD(loops)
# loop_refine.apply(pose)  # takes too long


# KIC
loops = protocols.loops.Loops()
loops.add_loop(loop1)

set_single_loop_fold_tree(pose, loop1)

sw_low = SwitchResidueTypeSetMover("centroid")
sw_low.apply(pose)
kic_perturb = protocols.loops.loop_mover.perturb.LoopMover_Perturb_KIC(loops)
# kic_perturb.apply(pose)  # won't show in Pymol for efficiency # takes too long


sw_high = SwitchResidueTypeSetMover("fa_standard")
sw_high.apply(pose)
kic_refine = protocols.loops.loop_mover.refine.LoopMover_Refine_KIC(loops)
# kic_refine.apply(pose)  # won't show in Pymol for efficiency # takes too long

# use the KinematicMover explicitly in centroid stage
sw_low.apply(pose)
kic_mover = KinematicMover()
kic_mover.set_pivots(16, 20, 24)
kic_mover.apply(pose)
