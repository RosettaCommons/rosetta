#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   best_beta_backbones.py
## @brief  Scan through beta configurations, ruling out bump check violations
## @author Andrew Watkins

def atoms_within_angle(pose, resi, ai, aj):
	neighbors = pose.residue_type(resi).bonded_neighbor(ai)
	for neighb in neighbors:
		if pose.residue_type(resi).atoms_are_bonded(neighb, aj): return true
	return False

def bump_check(pose, prop_vdw):
	for im1 in xrange(pose.size()):
		ii = im1 + 1
		for aim1 in xrange(pose.residue_type(ii).natoms());
			ai = aim1 + 1
			ai_r = pose.residue_type( ii ).atom_type( ai ).lj_radius();

			for jm1 in xrange(pose.size()):
				if jm1 < im1: continue
				jj = jm1 + 1
			
				for ajm1 in xrange(pose.residue_type(jj).natoms()):
					aj = ajm1 + 1
					if jj == ii and ai == aj: continue
					if jj == ii and pose.residue_type(ii).atoms_are_bonded(ai, aj): continue
					if jj == ii and atoms_within_angle(pose, ii, ai, aj): continue

					Real aj_r = pose.residue_type(jj).atom_type(aj).lj_radius();

					# skip bonded -- explicit for now!
					# Just make it so res 2's C can't clash with residue 3
					if ii == 3 and pose.residue(jj).atom_name(aj) == " C  ": continue
					if ii == 1 and pose.residue(jj).atom_name(aj) == " N  ": continue
					if jj == 3 and pose.residue(ii).atom_name(ai) == " C  ": continue
					if jj == 1 and pose.residue(ii).atom_name(ai) == " N  ": continue

					# Also handle 2 bond inter cases: 3-H 2-C, 3-N 2-O
					if ii == 1 and pose.residue_type( ii ).atom_name( ai ) == " C  " and jj == 2:
						neigh = pose.residue_type( jj ).bonded_neighbor( pose.residue_type( jj ).atom_index( " N  " ) )
						if aj in neigh: continue

					if ii == 2 and pose.residue_type( ii ).atom_name( ai ) == " N  " and jj == 1:
						neigh = pose.residue_type( jj ).bonded_neighbor( pose.residue_type( jj ).atom_index( " C  " ) )
						if aj in neigh: continue


					if ii == 2 and pose.residue_type( ii ).atom_name( ai ) == " C  " and jj == 3:
						neigh = pose.residue_type( jj ).bonded_neighbor( pose.residue_type( jj ).atom_index( " N  " ) )
						if aj in neigh: continue

					if ii == 3 and pose.residue_type( ii ).atom_name( ai ) == " N  " and jj == 2:
						neigh = pose.residue_type( jj ).bonded_neighbor( pose.residue_type( jj ).atom_index( " C  " ) )
						if aj in neigh: continue

					dsq = pose.residue(ii ).xyz(ai).distance_squared(pose.residue(jj).xyz(aj))
					if dsq < pow( ai_r + aj_r, 2 ) * prop_vdw * prop_vdw:
						return True

	return False

def report_proportion_accessible(rtname, bump_frac):
	print "Examining ", rtname, "..."

	score_fxn = core.scoring.get_score_function()
	pose = pyrosetta.pose_from_sequence("X[ACE]X[{}]X[NME]".format(rtname))

	# There are going to be 72^3 measurements.
	denominator = 72 * 72 * 72;
	numerator = 0;

	mm = kinematics.MoveMap()
	mm.set_chi(2, True);
	for ii in xrange(pose.residue_type(2).nchi()):
		pose.set_torsion(id.TorsionID(2, id.CHI, ii), 180);
	minm = protocols.minimization_packing.MinMover(mm, score_fxn, "linmin_iterated", 0.001, True)

	for phi in xrange(-175, 180, 5):
		for tht in xrange(-175, 180, 5):
			for psi in xrange(-175, 180, 5):
				pose.set_torsion(id.TorsionID(2, id.BB, 1), phi);
				pose.set_torsion(id.TorsionID(2, id.BB, 2), tht);
				pose.set_torsion(id.TorsionID(2, id.BB, 3), psi);
				min.apply(pose);
				if not bump_check(pose, bump_frac):
					print phi, tht, psi
					numerator += 1
	print "Proportion for", rtname, "is", numerator, "/", denominator, "or", float(numerator)/float(denominator)

def report_proportion_accessible_alpha(rtname, bump_frac):
	print "Examining ", rtname, "..."

	score_fxn = core.scoring.get_score_function()
	pose = pyrosetta.pose_from_sequence("X[ACE]X[{}]X[NME]".format(rtname))

	# There are going to be 72^3 measurements.
	denominator = 72 * 72;
	numerator = 0;

	mm = kinematics.MoveMap()
	mm.set_chi(2, True);
	minm = protocols.minimization_packing.MinMover(mm, score_fxn, "linmin_iterated", 0.001, True)

	for phi in xrange(-175, 180, 5):
		for psi in xrange(-175, 180, 5):
			pose.set_torsion(id.TorsionID(2, id.BB, 1), phi);
			pose.set_torsion(id.TorsionID(2, id.BB, 2), psi);
			min.apply(pose);
			if not bump_check(pose, bump_frac):
				print phi, psi
				numerator += 1
	print "Proportion for", rtname, "is", numerator, "/", denominator, "or", float(numerator)/float(denominator)


def report_proportions_accessible(bump_frac):

	beta_rts = ["B3A", "B3C", "B3D", "B3E", "B3F", "B3G", "B3H", "B3I", "B3K", "B3L", "B3M", "B3N", "B3P", "B3Q", "B3R", "B3S", "B3T", "B3V", "B3W", "B3Y"]

	for rt in beta_rts:
		report_proportion_accessible(rt, bump_frac)

	alpha_rts = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]

	for rt in alpha_rts:
		report_proportion_accessible_alpha(rt, bump_frac)

# Pass one argument: the bump-fraction of LJ radius below which you count a violation.
if __name__ == '__main__':
	report_proportions_accessible(sys.argv[1])
