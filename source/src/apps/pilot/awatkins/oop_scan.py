#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   oop_scan.py
## @brief  Creates an OOP dimer and toys around with its dihedrals.
## @author Andrew Watkins

def margin(val, comp, range):
	return val > comp - range and val < comp + range

if __name__ == '__main__':
	score_fxn = scoring.get_score_function()

	score_fxn.set_weight(core.scoring.atom_pair_constraint, 1.0)
	score_fxn.set_weight(core.scoring.angle_constraint, 1.0)
	score_fxn.set_weight(core.scoring.dihedral_constraint, 1.0)

	filename = sys.argv[1]
	pose = Pose()
	import_pose.pose_from_file( pose, filename , core.import_pose.PDB_file)

	anything_okay_for_orientation = False
	mrgn = 30

	while not anything_okay_for_orientation:

		for ( Size ii = 1 ii <= pose.size()-1 ++ii ) {
			if pose.residue_type(ii    ).name3() == "PRO": continue
			if pose.residue_type(ii + 1).name3() == "PRO": continue
			if pose.residue_type(ii    ).has_variant_type(OOP_POST): continue

			# extrapolate where the OOP carbon positions would be
			if pose.residue(ii).has("H") and pose.residue(ii + 1).has("H"):
				v1 = pose.residue(ii    ).atom("H").xyz() - pose.residue(ii    ).atom("N").xyz()
				v2 = pose.residue(ii + 1).atom("H").xyz() - pose.residue(ii + 1).atom("N").xyz()
				if (v2 - v1).length() < 1.6:
					print ii, "and", (ii + 1), "okay for distance"
					if margin( pose.phi(ii), -145, mrgn) and margin(pose.psi(ii), -10, mrgn)
							and (margin(pose.phi(ii + 1), -135, mrgn) or margin(pose.phi(ii + 1), 70, mrgn))
							and margin(pose.psi(ii + 1), 75, mrgn):
						print ii, " and ", (ii + 1), " okay for orientation"
						anything_okay_for_orientation = True

					if (margin(pose.phi(ii), 145, mrgn) and margin(pose.psi(ii), 10, mrgn)
							and (margin(pose.phi(ii + 1), 135, mrgn) or margin(pose.phi(ii + 1 , -70, mrgn))
							and margin(pose.psi(ii + 1), -75, mrgn)) {
						print ii, " and ", (ii + 1), " okay for orientation as Ds"
						anything_okay_for_orientation = True
			mrgn += 1

	print "Final required margin was", mrgn

