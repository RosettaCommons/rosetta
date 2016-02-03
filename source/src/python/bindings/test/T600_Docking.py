# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @author Sergey Lyskov

from rosetta import *
from rosetta.protocols.rigid import *
rosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

print 'Docking ----------------------------------------------------'

dock_p = pose_from_file("../test/data/test_dock.pdb")
dock_jump = 1
#DockingProtocol().setup_foldtree(dock_p)

to_centroid = protocols.simple_moves.SwitchResidueTypeSetMover('centroid')

jmp_arr = utility.vector1_int()
#jmp_arr = utility.vector1_ulong()
jmp_arr.append(1)
setup_foldtree(dock_p, '_', jmp_arr)

starting_p = Pose()
starting_p.assign(dock_p)

to_centroid.apply(dock_p)

dock_pert = RigidBodyPerturbMover(dock_jump, 3, 8)
dock_pert.apply(dock_p)

spin = RigidBodySpinMover( dock_jump )
spin.apply(dock_p)

slide_into_contact = DockingSlideIntoContact( dock_jump )
slide_into_contact.apply(dock_p)

docking_lowres = DockingLowRes()
docking_lowres.apply(dock_p)

#DockingProtocol().recover_sidechains(dock_p, starting_p)
recover_side_chain_mover = protocols.simple_moves.ReturnSidechainMover(starting_p)
recover_side_chain_mover.apply(dock_p)

#docking_highres = DockingHighRes()
#docking_highres.apply(dock_p)

# NEEDED: access and print interface
# TODO: Fix print fold_tree to avoid -1/1 codes
# TODO: Make jump_num default to 1 in all docking movers
