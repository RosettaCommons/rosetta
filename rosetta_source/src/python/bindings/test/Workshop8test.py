#! /usr/bin/python
# List of commands used in PyRosetts Workshop #2

from rosetta import *
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
score = create_score_function("standard")
jd = PyJobDistributor("output", 1, score)

lrms = loop_rmsd(pose, reference_pose, loops, True)
jd.additional_decoy_info = " LRMSD: " + str(lrms)

# High-Resolution Loop Protocol
loop_refine = LoopMover_Refine_CCD(loops)

# loop_refine.apply(pose)  # takes too long
