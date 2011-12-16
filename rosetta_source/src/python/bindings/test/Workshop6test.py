#! /usr/bin/python
# List of commands used in PyRosetts Workshop #2

from rosetta import *
init()

# Side Chain Conformations, the Rotamer Library, and Dunbrack Energies
pose = pose_from_pdb("test/data/workshops/1YY8.clean.pdb")
scorefxn = create_score_function("standard")

for i in range(1, 3):
    pose.set_chi(i, 49, 180)

# Monte Carlo Side-Chain Packing
task_pack = standard_packer_task(pose)
task_pack.restrict_to_repacking()
task_pack.temporarily_fix_everything()
task_pack.temporarily_set_pack_residue(49, True)

print task_pack

pack_mover = PackRotamersMover(scorefxn, task_pack)

pack_mover.apply(pose)

# Design
generate_resfile_from_pdb("test/data/workshops/1YY8.clean.pdb", "test/data/workshops/1YY8.resfile")
generate_resfile_from_pose(pose, "test/data/workshops/1YY8.resfile")

task_design = TaskFactory.create_packer_task(pose)
parse_resfile(pose, task_design, "test/data/workshops/1YY8.resfile")
