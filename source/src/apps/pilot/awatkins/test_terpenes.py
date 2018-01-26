#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   test_terpenes.py
## @brief  Test terpene RTs
## @author Andrew Watkins

if __name__ == '__main__':

	chm = rosetta.core.chemical.ChemicalManager.get_instance()
    rts = chm.residue_type_set('fa_standard').get_self_ptr()#.get_self_weak_ptr()
    
	score_fxn = get_score_function()

	Pose pose
	#make_pose_from_sequence(pose, "X[isoprene:terpene_lower_full]X[isoprene]X[isoprene]X[isoprene]X[isoprene]X[isoprene]X[isoprene]X[isoprene:terpene_upper_full]", rts)
	lower = conformation.ResidueFactory.create_residue(rts.name_map("isoprene:terpene_lower_full"))
	inner = conformation.ResidueFactory.create_residue(rts.name_map("isoprene"))
	upper = conformation.ResidueFactory.create_residue(rts.name_map("isoprene:terpene_upper_full"))

	pose.append_residue_by_jump(lower, 1, "", "", True)
	pose.append_residue_by_bond(inner, True)
	pose.append_residue_by_bond(inner, True)
	pose.append_residue_by_bond(inner, True)
	pose.append_residue_by_bond(inner, True)
	pose.append_residue_by_bond(inner, True)
	pose.append_residue_by_bond(upper, True)

	kinematics.MoveMapOP mm( new kinematics.MoveMap )
	mm.set_bb( true )
	mm.set_chi( true )

	pose.dump_pdb( "done.pdb" )
	for ii in xrange(pose.size()):
		i = ii + 1
		pose.set_torsion(id.TorsionID(i, id.BB, pose.residue(i).mainchain_torsions().size()-1), 60)
		pose.set_torsion(id.TorsionID(i, id.BB, pose.residue(i).mainchain_torsions().size()), 180)

	pose.dump_pdb( "set.pdb" )

	before_score = ( *score_fxn )( pose )
	#Real best_score = before_score
	mcpose = Pose(pose)
	naccept, nthermal, nreject = 0, 0, 0
	for i in xrange(100000):
		resi = int(numeric.random.rg().uniform() * pose.total_residue() + 1)
		# skip torsion 1, which is immobile.
		tori = int(numeric.random.rg().uniform() * (pose.residue(resi).mainchain_torsions().size() - 1) + 2)

		id.TorsionID torid( resi, id.BB, tori )
		Real torval = mcpose.torsion( torid )
		mcpose.set_torsion( torid, torval + numeric.random.rg().gaussian() * 10 )
		Real mcscore = ( *score_fxn )( mcpose )
		if mcscore < before_score:
			before_score = mcscore
			#best_score = mcscore
			pose = mcpose
			#TR, "Accepted.", std.endl
			++naccept
		elif float(exp( -( mcscore - before_score ) )) > numeric.random.rg().uniform():
			print "Accepted thermally (scores ", before_score, " then ", mcscore, " prob ", Real(exp( -( mcscore - before_score ) ) ), ")."
			before_score = mcscore
			pose = Pose(mcpose)
			nthermal += 1
		else:
			#TR, "Rejected.", std.endl
			nreject += 1
		
		if i % 1000 == 0:
			TR, "Report.", std.endl
			TR, i, " moves. ", naccept, " accepted ( ", ( 100 * Real( naccept/i ) ), " ).", std.endl
			TR, nthermal, " accepted thermally ( ", ( 100 * Real( nthermal/i ) ), " ).  "
			TR, nreject, " rejected ( ", ( 100 * Real( nreject/i ) ), " ).", std.endl


	pose.dump_pdb("mc.pdb")

	min = MinMover(mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, True)
	min.apply(pose)

	pose.dump_pdb("min.pdb")
