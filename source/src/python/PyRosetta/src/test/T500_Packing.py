# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @author Sergey Lyskov

from __future__ import print_function

from rosetta import *
from pyrosetta import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


print('Packing and Design ----------------------------------------------')

print('mover: PackRotamersMover')
pose = core.import_pose.pose_from_file("../test/data/test_in.pdb")
scorefxn = get_fa_scorefxn() #  create_score_function('standard')
#task_pack = core.pack.task.TaskFactory.create_packer_task(pose)
task_pack = standard_packer_task(pose)
task_pack.restrict_to_repacking()
task_pack.nonconst_residue_task(5).prevent_repacking()
print(task_pack)
pack = protocols.simple_moves.PackRotamersMover( scorefxn, task_pack )
pack.apply(pose)

rotamer_trials = protocols.simple_moves.RotamerTrialsMover(scorefxn, task_pack)
rotamer_trials.apply(pose)

#TODO add task interface commands like:
#prevent_repacking(1,total_residue) or (res) or (True [=all])
#allow_repacking() # same argument options/overload
# ... needs discussion with Mini community...task is 'commutative'

#TODO add meaningful info to print pack
print(pack)

task_design = core.pack.task.TaskFactory.create_packer_task(pose)

# mjo -> in refactoring the resfile reader I modified this to use parse_resfile
#task_design.read_resfile("../test/data/test_in.resfile")
core.pack.task.parse_resfile(pose, task_design, "../test/data/test_in.resfile")# TODO this hard-crashes if file is not a resfile


# BUG reported on jenkins machine: comments before 'start' line cause crash
print(task_design)
pack = protocols.simple_moves.PackRotamersMover( scorefxn, task_design )
pack.apply(pose)

#TaskFactory options
tf = standard_task_factory()
tf.push_back(core.pack.task.operation.RestrictToRepacking())

pr = core.pack.task.operation.PreventRepacking()
pr.include_residue(5)

tf.push_back(pr)

pack = protocols.simple_moves.PackRotamersMover( scorefxn )
pack.task_factory(tf)
pack.apply(pose)
