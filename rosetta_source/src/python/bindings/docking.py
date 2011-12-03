# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
from rosetta import *
init()

p = Pose()
pose_from_pdb(p, "test_dock.pdb")

starting_p = Pose()
starting_p.assign(p) #saving the starting structure

print "setting up docking fold tree"
dock_prot = DockingProtocol() #DockingProtocol object contains many useful docking functions
dock_prot.setup_foldtree(p)
dock_jump = 1

print "set up scoring functions"
scorefxn_low = create_score_function('interchain_cen')
scorefxn_high = create_score_function('docking')
scorefxn_high_min = create_score_function_ws_patch('docking','docking_min')

print "setting up movers"

#centroid/fullatom conversion movers
to_centroid = protocols::simple_moves::SwitchResidueTypeSetMover('centroid')
to_fullatom = protocols::simple_moves::SwitchResidueTypeSetMover('fa_standard')
recover_sidechains = protocols::simple_moves::ReturnSidechainMover(starting_p)

#initial perturbation movers
randomize1 = RigidBodyRandomizeMover(p, dock_jump, rigid::partner_upstream)
randomize2 = RigidBodyRandomizeMover(p, dock_jump, rigid::partner_downstream)
dock_pert = RigidBodyPerturbMover(dock_jump, 3, 8) #3A translation, 8 degrees rotation
spin = RigidBodySpinMover( dock_jump )
slide_into_contact = DockingSlideIntoContact( dock_jump )

#docking lowres movers
docking_lowres = DockingLowRes( scorefxn_low, dock_jump )

#docking highres movers
docking_highres = DockingHighRes( scorefxn_high_min, dock_jump )

print "set up job distributor"
jd = PyJobDistributor("dock_output", 20, scorefxn_high)
jd.native_pose = starting_p

print "beginning docking..."

print "convert to centroid mode"
to_centroid.apply(p)
starting_p_centroid = Pose()
starting_p_centroid.assign(p)

while (jd.job_complete == False):
	p.assign(starting_p_centroid)

	print "initial perturbation"
	dock_pert.apply(p)
	slide_into_contact.apply(p)

	print "low resolution stage docking"
	docking_lowres.apply(p)

	print "convert to fullatom mode"
	to_fullatom.apply(p)
#	recover_sidechains.apply(p)
	#dock_prot.recover_sidechains(p)  # <-- requare native pose, C++ function changed signature, commenting out for now

	print "high resolution stage docking"
	docking_highres.apply(p)

	print "outputting decoy..."
	jd.output_decoy(p)

print "docking complete!"
