#! /usr/bin/python
# List of commands used in PyRosetts Workshop #7

import sys

from rosetta import *
import rosetta.protocols.rigid as rigid_moves

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


'''
for _i in range(10):
    try:
'''

# Docking Moves in Rosetta
pose = pose_from_file("../test/data/workshops/complex.start.pdb")

print pose.fold_tree()

setup_foldtree(pose, "A_B", Vector1([1]))
print pose.fold_tree()

jump_num = 1
print pose.jump(jump_num).get_rotation()
print pose.jump(jump_num).get_translation()

print "_____ Check point 1"
pert_mover = rigid_moves.RigidBodyPerturbMover(jump_num, 8, 3)
#pert_mover.apply(pose)

randomize1 = rigid_moves.RigidBodyRandomizeMover(pose, jump_num, rigid_moves.partner_upstream)
randomize2 = rigid_moves.RigidBodyRandomizeMover(pose, jump_num, rigid_moves.partner_downstream)

print "_____ Check point 2"
#randomize1.apply(pose)
#randomize2.apply(pose)
slid = DockingSlideIntoContact(jump_num)
slide = FaDockingSlideIntoContact(jump_num)
slide.apply(pose)

movemap = MoveMap()
movemap.set_jump(jump_num, True)

scorefxn = create_score_function("talaris2013")
scorefxn( pose )

print "_____ Check point 3"
print 'Making MinMover...'
min_mover = MinMover()
min_mover.movemap(movemap)
min_mover.score_function(scorefxn)

#min_mover.apply(pose)

print 'Done Applying MinMover!'


'''
    #except rosetta.PyRosettaException: pass
    except RuntimeError: pass

else:
    print 'Was not able to finish min_mover in 10 tries, failing...'
    sys.exit(1)
'''

# Low-Resolution Docking via RosettaDock
switch_low = SwitchResidueTypeSetMover("centroid")
pose_high = Pose()
pose_high.assign(pose)

switch_low.apply(pose)
pose_low = Pose()
pose_low.assign(pose)

setup_foldtree(pose_low, "A_B", Vector1([1]))

scorefxn_low = create_score_function("interchain_cen")

dock_lowres = DockingLowRes(scorefxn_low, jump_num)
dock_lowres.apply(pose_low)

print CA_rmsd(pose, pose_low)
print calc_Lrmsd(pose, pose_low, Vector1([1]))

# Job Distributor
import tempfile
output = tempfile.mkstemp()[1]

jd = PyJobDistributor(output, 10, scorefxn_low)

native_pose = pose_from_file("../test/data/workshops/complex.high.pdb")
jd.native_pose = native_pose

starting_pose = Pose()
starting_pose.assign(pose_low)

while (jd.job_complete == False):
    pose_low.assign(starting_pose)
    dock_lowres.apply(pose_low)
    jd.output_decoy(pose_low)

# High-Resolution Docking
scorefxn_high = create_score_function_ws_patch("pre_talaris_2013_standard.wts", "docking")
dock_hires = DockMCMProtocol()
dock_hires.set_scorefxn(scorefxn_high)
dock_hires.set_partners("A_B")

recover_sidechains = ReturnSidechainMover(pose_high)
recover_sidechains.apply(pose)

print "done"
