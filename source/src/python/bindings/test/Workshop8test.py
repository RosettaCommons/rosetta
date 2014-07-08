#! /usr/bin/python
# List of commands used in PyRosetts Workshop #2

from rosetta import *
from rosetta.protocols.loops.loop_closure.ccd import *
from rosetta.protocols.loops.loop_mover.refine import *
from rosetta.protocols.loops.loop_mover.perturb import *
from rosetta.protocols.loops.loop_closure.kinematic_closure import *

init()

pose = pose_from_pdb("test/data/test_in.pdb")

# Fold Tree
ft = FoldTree()
ft.add_edge(1, 13, -1)
ft.add_edge(13, 19, -1)
ft.add_edge(13, 26, 1)
ft.add_edge(26, 20, -1)
ft.add_edge(26, 116, -1)

print ft
ft.check_fold_tree()

pose.fold_tree(ft)

for res in (10, 13, 16, 23, 26, 30):
    pose.set_phi(res, 180)
    pose.dump_pdb("loop" + str(res) + ".pdb")

pmm = PyMOL_Mover()
pmm.apply(pose)
pmm.send_foldtree(pose)
pmm.view_foldtree_diagram(pose, ft)

ft.clear()
ft.simple_tree(116)
ft.new_jump(76, 85, 80)

# Cyclic Coordination Descent (CCD) Loop Closure
movemap = MoveMap()
movemap.set_bb(True)
movemap.set_chi(True)

loop1 = Loop(15, 24, 19)
ccd = CcdLoopClosureMover(loop1, movemap)

ccd.apply(pose)
set_single_loop_fold_tree(pose, loop1)

# Multiple Loops
loop2 = Loop(78, 83, 80)

loops = Loops()
loops.add_loop(loop1)
loops.add_loop(loop2)

# Loop Building
reference_pose = pose_from_pdb("test/data/test_in.pdb")
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
loops = Loops()
loops.add_loop(loop1)

set_single_loop_fold_tree(pose, loop1)

sw_low = SwitchResidueTypeSetMover("centroid")
sw_low.apply(pose)
kic_perturb = LoopMover_Perturb_KIC(loops)
# kic_perturb.apply(pose)  # won't show in Pymol for efficiency # takes too long


sw_high = SwitchResidueTypeSetMover("fa_standard")
sw_high.apply(pose)
kic_refine = LoopMover_Refine_KIC(loops)
# kic_refine.apply(pose)  # won't show in Pymol for efficiency # takes too long

# use the KinematicMover explicitly in centroid stage
sw_low.apply(pose)
kic_mover = KinematicMover()
kic_mover.set_pivots(16, 20, 24)
kic_mover.apply(pose)




