#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   pna.py
## @brief  Look for structure in peptide nucleic acids
## @author Andrew Watkins

def add_hbond_constraints(pose):

	HarmonicFuncOP harm( new core.scoring.func.HarmonicFunc( 2, .5 ) )
	CircularHarmonicFuncOP pi( new core.scoring.func.CircularHarmonicFunc( 3.14159, .05 ) )
	CircularHarmonicFuncOP zero( new core.scoring.func.CircularHarmonicFunc( 0, .05 ) )

	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
		*new AtomID( pose.residue( 1 ).atom_index( "NL1" ), 1 ),
		*new AtomID( pose.residue( 4 ).atom_index( "1HI2" ), 4 ), harm ) ) )
	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
		*new AtomID( pose.residue( 1 ).atom_index( "1HL2" ), 1 ),
		*new AtomID( pose.residue( 4 ).atom_index( "OL" ), 4 ) , harm ) ) )

	pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
		*new AtomID( pose.residue( 4 ).atom_index( "OL" ), 4 ),
		*new AtomID( pose.residue( 1 ).atom_index( "1HL2" ), 1 ),
		*new AtomID( pose.residue( 1 ).atom_index( "2HL2" ), 1 ),
		*new AtomID( pose.residue( 1 ).atom_index( "CK2" ), 1 ), pi ) ) )

	pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
		*new AtomID( pose.residue( 4 ).atom_index( "1HI2" ), 4 ),
		*new AtomID( pose.residue( 1 ).atom_index( "NL1" ), 1 ),
		*new AtomID( pose.residue( 1 ).atom_index( "CK1" ), 1 ),
		*new AtomID( pose.residue( 1 ).atom_index( "CK2" ), 1 ), pi ) ) )
	
if __name__ == '__main__':
		
	pose = Pose()
	scorefxn = get_score_function()
	guaranteed_cart_scorefxn = ScoreFunctionFactory.create_score_function("talaris2013_cart")
	#scorefxn.set_weight( hbond_sc, 5 )

	#Get the residue set we are drawing from.
	chm = rosetta.core.chemical.ChemicalManager.get_instance()
    rts = chm.residue_type_set('fa_standard').get_self_ptr()#.get_self_weak_ptr()

	pna = residue_set_cap.name_map("APN")
	pnc = residue_set_cap.name_map("CPN")
	png = residue_set_cap.name_map("GPN")
	pnt = residue_set_cap.name_map("TPN")
	pnu = residue_set_cap.name_map("UPN")

	res_pna_nterm = Residue(residue_set_cap.name_map("APN:NtermProteinFull"), True)
	res_pna_cterm = Residue(residue_set_cap.name_map("APN:CtermProteinFull"), True)
	res_pna = Residue(pna, True)
	res_pnc = Residue(pnc, True)
	res_png = Residue(png, True)
	res_pnt = Residue(pnt, True)
	res_pnu = Residue(pnu, True)
	res_pnu_nterm = Residue(residue_set_cap.name_map("UPN:NtermProteinFull"), True)
	res_pnu_cterm = Residue(residue_set_cap.name_map("UPN:CtermProteinFull"), True)

	pose.append_residue_by_jump(res_pna, 1)
	pose.append_residue_by_bond(res_pnc, True)
	pose.append_residue_by_jump(res_png, 2)
	pose.append_residue_by_bond(res_pnu, True)

	trans_mover = protocols.rigid.RigidBodyTransMover(pose, 1)
	trans_mover.step_size(4)
	trans_mover.apply(pose)
	big_pert_dock_rbpm = rigid.RigidBodyPerturbMover(1, 1.0, 1.5)
	pert_dock_rbpm = rigid.RigidBodyPerturbMover(1, 1.0, 0.5)

	add_hbond_constraints(pose)

	mm = kinematics.MoveMap()
	mm.set_bb(True)

	min_mm = kinematics.MoveMap()
	#min_mm.set_chi(True)
	min_mm.set_jump(True)

	# SET INITIAL CONFORMATION
	for resim1 in xrange(pose.size()):
		resi = resim1 + 1
		mm.set(id.TorsionID(resi, id.BB, 5), False)
		min_mm.set(id.TorsionID(resi, id.BB, 5), False)

		pose.set_torsion(id.TorsionID(resi, id.BB, 1), 70)
		pose.set_torsion(id.TorsionID(resi, id.BB, 2), 70)
		pose.set_torsion(id.TorsionID(resi, id.BB, 3), 95)
		pose.set_torsion(id.TorsionID(resi, id.BB, 4), -5)
		pose.set_torsion(id.TorsionID(resi, id.BB, 5), 180)

		pose.set_torsion(id.TorsionID(resi, id.CHI, 1), 0)
		pose.set_torsion(id.TorsionID(resi, id.CHI, 2), 180)
		pose.set_torsion(id.TorsionID(resi, id.CHI, 3), -90)

	small_tor_mover = RandomTorsionMover(mm, 1, 10)
	min = MinMover(min_mm, guaranteed_cart_scorefxn, "lbfgs_armijo_nonmonotone", 0.001, True)
	min.cartesian(True)
	pert_mc = moves.MonteCarlo(pose, scorefxn, 3)

	min.apply(pose)
	pose.dump_pdb("min.pdb")

	apc = 0.1
	scorefxn.set_weight(atom_pair_constraint, apc)
	scorefxn.set_weight(angle_constraint, apc)
	scorefxn.set_weight(dihedral_constraint, apc)

	best_score = 100000
	for ii in xrange(1000):
		print "PNA Creator: Round", ii, "/ 1000"
		print "PNA Creator: score", scorefxn(pose), "best", best_score

		for jj in xrange(10):
			small_tor_mover.apply( pose )
			
		if pert_mc.boltzmann(pose):
			min.apply(pose)
			test_score = pose.energies().total_energies().dot(scorefxn.weights())
			if test_score <= best_score:
				best_score = test_score
				pert_mc.reset(pose)

		if ii % 100 == 0: 
			pose.dump_pdb("out_{}.pdb".format(ii))
			apc += 0.1


	print "about to dump"
	pose.dump_pdb("out.pdb")

	print "PNA Creator: Removing constraints and remodeling"
	pose.constraint_set().clear()

	best_score = 100000
	for ii in xrange(100):
		print "PNA Creator: Round", ii, "/ 100"
		print "PNA Creator: score", scorefxn(pose), "best", best_score

		for jj in xrange(100):
			small_tor_mover.apply( pose )
			
		if pert_mc.boltzmann(pose):
			test_score = pose.energies().total_energies().dot(scorefxn.weights())
			if test_score <= best_score:
				min.apply( pose )
				best_score = test_score
				pert_mc.reset(pose)

		if ii % 10 == 0:
			pose.dump_pdb("new_{}.pdb".format(ii))

	pert_mc.recover_low(pose)
	pose.dump_pdb("new.pdb")
