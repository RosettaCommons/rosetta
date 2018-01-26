#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   test_nu_angles.py
## @brief  mmm
## @author Andrew Watkins

if __name__ == '__main__':

	#devel::init(argc, argv)

	score_fxn = get_score_function()
	score_fxn.set_weight(ring_close, 1)

	pose = pyrosetta.pose_from_sequence("X[A04]X[A98]X[B02]X[B06]X[C01]X[B19]X[B95]X[C00]X[C12]")

	for iim1 in xrange(pose.size()):
		ii = iim1 + 1
		pose.set_phi(  ii, -150)
		pose.set_psi(  ii,  150)
		pose.set_omega(ii,  180)
	
	pose.dump_scored_pdb("out.pdb", score_fxn)

	mm = kinematics.MoveMap()
	mm.set_bb(True)
	mm.set_chi(True)
	mm.set_nu(True)

	min = MinMover(mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.0001, True)
	min.apply(pose)

	pose.dump_scored_pdb("min.pdb", score_fxn)