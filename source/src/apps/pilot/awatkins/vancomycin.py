#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   vancomycin.py
## @brief  n/a
## @author Andrew Watkins

def add_constraints(vancomycin):
	CZ2_2(vancomycin.residue(2).type().atom_index("CZ2"), 2)
	OT_2(vancomycin.residue(2).type().atom_index("OT"), 2)
	CZ2_6(vancomycin.residue(6).type().atom_index("CZ2"), 6)
	OT_6(vancomycin.residue(6).type().atom_index("OT"), 6)
	CD1(vancomycin.residue(4).type().atom_index("CD1"), 4)
	CE(vancomycin.residue(4).type().atom_index("CE"), 4)
	CG1(vancomycin.residue(4).type().atom_index("CG1"), 4)
	CG2(vancomycin.residue(4).type().atom_index("CG2"), 4)
	CD2(vancomycin.residue(4).type().atom_index("CD2"), 4)
	CD12(vancomycin.residue(5).type().atom_index("CD1"), 5)
	CG12(vancomycin.residue(7).type().atom_index("CG1"), 7)
	CB(vancomycin.residue(7).type().atom_index("CB"), 7)
	CD13(vancomycin.residue(7).type().atom_index("CD1"), 7)

	harm = core.scoring.func.HarmonicFunc(1.31, 0.02)
	vancomycin.add_constraint(core.scoring.constraints.AtomPairConstraint(OT_2, CD1, harm))
	vancomycin.add_constraint(core.scoring.constraints.AtomPairConstraint(OT_6, CD2, harm))
	vancomycin.add_constraint(core.scoring.constraints.AtomPairConstraint(CD12, CG12, harm))

	aharm(core.scoring.func.CircularHarmonicFunc( 104.5*3.14159/180, 0.02)
	aharm2(core.scoring.func.CircularHarmonicFunc( 120  *3.14159/180, 0.02)

	vancomycin.add_constraint(core.scoring.constraints.AngleConstraint(CZ2_2, OT_2, CD1, aharm2))
	vancomycin.add_constraint(core.scoring.constraints.AngleConstraint(CZ2_6, OT_6, CD2, aharm2))
	vancomycin.add_constraint(core.scoring.constraints.AngleConstraint(OT_2, CD1, CE, aharm2))
	vancomycin.add_constraint(core.scoring.constraints.AngleConstraint(OT_6, CD2, CE, aharm2))
	vancomycin.add_constraint(core.scoring.constraints.AngleConstraint(CD12, CG12, CD13, aharm2))

	charm = core.scoring.func.CircularHarmonicFunc(3.14, 0.02))
	vancomycin.add_constraint(core.scoring.constraints.DihedralConstraint(OT_2, CD1, CG1, CE, charm))
	vancomycin.add_constraint(core.scoring.constraints.DihedralConstraint(OT_6, CD2, CG2, CE, charm))
	vancomycin.add_constraint(core.scoring.constraints.DihedralConstraint(CD12, CG12, CB, CD13, charm))

	dharm = core.scoring.func.CircularHarmonicFunc(-71.5 * 3.14159 / 180, 0.02))
	OZ4 = core.id.AtomID(vancomycin.residue(4).type().atom_index("OZ"), 4)
	C18 = core.id.AtomID(vancomycin.residue(8).type().atom_index("C1"), 8)
	C28 = core.id.AtomID(vancomycin.residue(8).type().atom_index("C2"), 8)
	O28 = core.id.AtomID(vancomycin.residue(8).type().atom_index("O2"), 8)
	C19 = core.id.AtomID(vancomycin.residue(9).type().atom_index("C1"), 9)
	O58 = core.id.AtomID(vancomycin.residue(8).type().atom_index("O5"), 8)
	O59 = core.id.AtomID(vancomycin.residue(9).type().atom_index("O5"), 9)
	bond4 = core.scoring.constraints.AtomPairConstraint(OZ4, C18, harm)
	bond5 = core.scoring.constraints.AtomPairConstraint(O28, C19, harm)
	ang6 = core.scoring.constraints.AngleConstraint(CE, OZ4, C18, aharm2)
	ang7 = core.scoring.constraints.AngleConstraint(C28, O28, C19, aharm2)

	dih4 = core.scoring.constraints.DihedralConstraint(CE, OZ4, C18, O58, dharm)
	dih5 = core.scoring.constraints.DihedralConstraint(C28, O28, C19, O59, dharm)

	vancomycin.add_constraint(bond4)
	vancomycin.add_constraint(bond5)
	vancomycin.add_constraint(ang6)
	vancomycin.add_constraint(ang7)

	vancomycin.add_constraint(dih4)
	vancomycin.add_constraint(dih5)


if __name__ == '__main__':
	#option[ chemical.patch_selectors ].push_back( "CTERM_AMIDATION" )

	devel.init(argc, argv)

	scorefxn = get_score_function()
	scorefxn.set_weight_if_zero(atom_pair_constraint, 0.1)
	scorefxn.set_weight_if_zero(angle_constraint, 1.0)
	scorefxn.set_weight_if_zero(dihedral_constraint, 1.0)
	
	# Vanc sequence: N-methyl D-leu
	chm = rosetta.core.chemical.ChemicalManager.get_instance()
    restype_set_cap = chm.residue_type_set('fa_standard').get_self_ptr()#.get_self_weak_ptr()
    
	dleu_type = residue_set_cap.name_map("DLEU:NtermProteinMethylated")
	#ResidueType const & dleu_type = residue_set_cap.name_map( "DTRP:NtermProteinMethylated" )
	V01_type = residue_set_cap.name_map("DV01:aryl-O-conjugated")
	asn_type = residue_set_cap.name_map("ASN")
	#ResidueType const & asn_type = residue_set_cap.name_map( "LYS" )
	V02_type = residue_set_cap.name_map("DV02:aryl-O-conjugated:phg_cd1_conjugation:phg_cd2_conjugation")
	V02_type2 = residue_set_cap.name_map("DV02:phg_cd1_conjugation")
	V04_type = residue_set_cap.name_map("V04:aryl-O-conjugated")
	V03_type = residue_set_cap.name_map("V03:CtermProteinFull:aryl-C-conjugated")
	#ResidueType const & V03_type = residue_set_cap.name_map( "V03:aryl-C-conjugated" )

	dleu = Residue(dleu_type, True)
	V01 = Residue(V01_type, True)
	asn = Residue(asn_type, True)
	V02 = Residue(V02_type, True)
	V022 = Residue(V02_type2, True)
	V04 = Residue(V04_type, True)
	V03 = Residue(V03_type, True)

	vancomycin = Pose()
	vancomycin.append_residue_by_jump(dleu, 1)
	vancomycin.append_residue_by_bond(V01, True)
	vancomycin.append_residue_by_bond(asn, True)
	vancomycin.append_residue_by_bond(V02, True)
	vancomycin.append_residue_by_bond(V022, True)
	vancomycin.append_residue_by_bond(V04, True)
	vancomycin.append_residue_by_bond(V03, True)
	#vancomycin.append_residue_by_jump( Residue( residue_set_cap.name_map( ".4)-beta-D-Glcp" ), true ), 2 )
	print "Appending by jump"
	vancomycin.append_residue_by_jump(Residue(residue_set_cap.name_map("->2)-beta-D-Glcp"), True), 2)
	#vancomycin.append_residue_by_bond( Residue( residue_set_cap.name_map( ".4)-Vnc:non-reducing_end:3-NH3+" ), true ), true )
	vancomycin.append_residue_by_bond(Residue(residue_set_cap.name_map("->4)-alpha-Daup3Me:non-reducing_end:3-Me"), True), True)

	translate = protocols.rigid.RigidBodyTransMover(vancomycin, 1)
	translate.step_size(4)
	translate.apply(vancomycin)

	# extra!
	#vancomycin.append_residue_by_bond( Residue( residue_set_cap.name_map( "GLN:CtermProteinFull" ), true ), true )
	#import_pose.set_reasonable_fold_tree( vancomycin )
	#FOLD_TREE  EDGE 1 2 -1  EDGE 2 7 -1  EDGE 2 8 1
	#FOLD_TREE  EDGE 1 2 -1  EDGE 2 7 -1  EDGE 4 8 -2

	#Size v02_connid1 = V02_type.residue_connection_id_for_atom( pose.residue( 4 ).atom_index( "CD1" ) )
	#Size v02_connid2 = V02_type.residue_connection_id_for_atom( pose.residue( 4 ).atom_index( "CD2" ) )
	#new_cys.residue_connection_partner( cys_connid, resi_vdp, vdp_connid )

	vancomycin.conformation().declare_chemical_bond(2, "OT", 4, "CD1")
	vancomycin.conformation().declare_chemical_bond(6, "OT", 4, "CD2")
	vancomycin.conformation().declare_chemical_bond(5, "CD1", 7, "CG1")
	vancomycin.conformation().declare_chemical_bond(4, "OZ", 8, "C1")
	vancomycin.conformation().declare_chemical_bond(8, "O2", 9, "C1")

	print "Declared bond"

	add_constraints(vancomycin)

	for i in xrange(vancomycin.size()-1):
		ii = i + 1
		vancomycin.set_phi(ii, -150)
		vancomycin.set_psi(ii, 150)
		vancomycin.set_omega(ii, 180)

	vancomycin.dump_pdb("vancomycin.pdb")

	pert_mm = kinematics.MoveMap()
	for i in xrange(vancomycin.size()):
		pert_mm.set_bb(i + 1, True)
		pert_mm.set_chi(i + 1, True)

	pert_mm.set_jump(1, True)
	#pert_mm.set_branches( 4, true )
	min_mover = minimization_packing.MinMover(pert_mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, True)

	for apc in xrange(0, 2, 0.01):
		scorefxn.set_weight(atom_pair_constraint, apc)
		scorefxn.set_weight(angle_constraint, apc)
		scorefxn.set_weight(dihedral_constraint, apc)

		min_mover.apply(vancomycin)

	vancomycin.dump_pdb("vancomycin_min1.pdb")

	f = vancomycin.fold_tree()
	f.clear()
	f.add_edge(1, 2, -1)
	f.add_edge(2, 7, -1)
	f.add_edge(4, 8, -2)
	f.add_edge(8, 9, -1)
	vancomycin.fold_tree(f)
	print "Set new fold tree ", vancomycin.fold_tree()
	pert_mm.set_branches(4, True)

	for apc in xrange(2, 0, -0.01):
		scorefxn.set_weight(atom_pair_constraint, apc)
		scorefxn.set_weight(angle_constraint, apc)
		scorefxn.set_weight(dihedral_constraint, apc)

		min_mover.apply(vancomycin)

	vancomycin.dump_pdb("vancomycin_min2.pdb")

	task = TaskFactory.create_packer_task(vancomycin)
	for ii in xrange(vancomycin.size()):
		i = ii + 1
		task.nonconst_residue_task(i).restrict_to_repacking()
		task.nonconst_residue_task(i).initialize_from_command_line()
	
	pack = PackRotamersMover(scorefxn, task)
	pack.apply(vancomycin)
	vancomycin.dump_pdb("vancomycin_pack.pdb")


	print "Just those branch res "
	vanc2 = Pose()
	vanc2.append_residue_by_jump(Residue(residue_set_cap.name_map("GLY:NtermProteinFull"), True), 1)
	vanc2.append_residue_by_bond(Residue(residue_set_cap.name_map("DV02:aryl-O-conjugated"), True), True)
	vanc2.append_residue_by_bond(Residue(residue_set_cap.name_map("GLY:CtermProteinFull"), True), True)
	print "Appending by jump"
	vanc2.append_residue_by_jump(Residue(residue_set_cap.name_map("->2)-beta-D-Glcp"), True), 2)
	vanc2.append_residue_by_bond(Residue(residue_set_cap.name_map("->4)-alpha-Daup3Me:non-reducing_end:3-Me"), True), True)

	vanc2.conformation().declare_chemical_bond(2, "OZ", 4, "C1")
	vanc2.conformation().declare_chemical_bond(4, "O2", 5, "C1")

	trans2 = protocols.rigid.RigidBodyTransMover(vanc2, 1)
	trans2.step_size(4)
	trans2.apply(vanc2)
	
	dharm = core.scoring.func.CircularHarmonicFunc(-71.5 * 3.14159 / 180, 0.02)
	dharm2 = core.scoring.func.CircularHarmonicFunc(176.9 * 3.14159 / 180, 0.02)
	dharm3 = core.scoring.func.CircularHarmonicFunc(60 * 3.14159 / 180, 0.02)

	harm = core.scoring.func.HarmonicFunc(1.31, 0.02)
	aharm2 = core.scoring.func.CircularHarmonicFunc(120 * 3.14159 / 180, 0.02)

	OZ4 = core.id.AtomID(vanc2.residue(2).type().atom_index("OZ"), 2)
	C18 = core.id.AtomID(vanc2.residue(4).type().atom_index("C1"), 4)
	C28 = core.id.AtomID(vanc2.residue(4).type().atom_index("C2"), 4)
	O28 = core.id.AtomID(vanc2.residue(4).type().atom_index("O2"), 4)
	C19 = core.id.AtomID(vanc2.residue(5).type().atom_index("C1"), 5)
	O58 = core.id.AtomID(vanc2.residue(4).type().atom_index("O5"), 4)
	O59 = core.id.AtomID(vanc2.residue(5).type().atom_index("O5"), 5)
	C58 = core.id.AtomID(vanc2.residue(4).type().atom_index("C5"), 4)
	C59 = core.id.AtomID(vanc2.residue(5).type().atom_index("C5"), 5)
	CE  = core.id.AtomID(vanc2.residue(2).type().atom_index("CE"), 2)

	vanc2.add_constraint(core.scoring.constraints.AtomPairConstraint(OZ4, C18, harm))
	vanc2.add_constraint(core.scoring.constraints.AtomPairConstraint(O28, C19, harm))
	vanc2.add_constraint(core.scoring.constraints.AngleConstraint(CE, OZ4, C18, aharm2))
	vanc2.add_constraint(core.scoring.constraints.AngleConstraint(C28, O28, C19, aharm2))

	vanc2.add_constraint(core.scoring.constraints.DihedralConstraint(CE, OZ4, C18, O58, dharm3))
	vanc2.add_constraint(core.scoring.constraints.DihedralConstraint(C28, O28, C19, O59, dharm))
	vanc2.add_constraint(core.scoring.constraints.DihedralConstraint(OZ4, C18, O58, C58, dharm2))

	pert_mm2 = kinematics.MoveMap()
	pert_mm2.set_bb(1, True)
	pert_mm2.set_chi(1, True)
	pert_mm2.set_bb(3, True)
	pert_mm2.set_chi(3, True)
	pert_mm2.set_jump(1, True)
	min_mover2 = minimization_packing.MinMover(pert_mm2, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, True)

	for apc in xrange(0, 2, 0.01):
		scorefxn.set_weight(atom_pair_constraint, apc)
		scorefxn.set_weight(angle_constraint, apc)
		scorefxn.set_weight(dihedral_constraint, apc)

		min_mover2.apply(vanc2)

	vanc2.dump_pdb("vanc2_min1.pdb")

	f2 = vanc2.fold_tree()
	f2.clear()
	f2.add_edge(1, 3, -1)
	f2.add_edge(2, 4, -2)
	f2.add_edge(4, 5, -1)
	vanc2.fold_tree(f2)
	print "Set new fold tree ", vanc2.fold_tree()
	pert_mm2.set_branches(2, True)
	pert_mm2.set_bb(2, True)
	pert_mm2.set_chi(2, True)
	pert_mm2.set_bb(4, True)
	pert_mm2.set_chi(4, True)
	pert_mm2.set_bb(5, True)
	pert_mm2.set_chi(5, True)

	for apc in xrange(2, 0, -0.01):
		scorefxn.set_weight(atom_pair_constraint, apc)
		scorefxn.set_weight(angle_constraint, apc)
		scorefxn.set_weight(dihedral_constraint, apc)

		min_mover2.apply(vanc2)

	vanc2.dump_pdb("vanc2_min2.pdb")

	task2 = TaskFactory.create_packer_task(vanc2)
	for i in xrange(vanc2.size()):
		task2.nonconst_residue_task(i + 1).restrict_to_repacking()
		task2.nonconst_residue_task(i + 1).initialize_from_command_line()

	pack2 = PackRotamersMover(scorefxn, task2)
	pack2.apply(vanc2)
	vanc2.dump_pdb("vanc2_pack.pdb")
