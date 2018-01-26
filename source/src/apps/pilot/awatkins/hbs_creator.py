#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   hbs_creator.py
## @brief  Creates an HBS from a peptide
## @author Andrew Watkins

# Gonna have to implement some argparse for options


# application specific options
#namespace hbs_creator {
# pert options
#StringOptionKey const hbs_chain ( "hbs_creator.hbs_chain" )
#IntegerOptionKey const hbs_final_res ( "hbs_creator.hbs_final_res" )
#IntegerOptionKey const hbs_length ( "hbs_creator.hbs_length" )
#BooleanOptionKey const final_repack( "hbs_creator.final_repack" )
#BooleanOptionKey const final_minimize( "hbs_creator.final_minimize" )
#BooleanOptionKey const final_mc ( "hbs_creator.final_mc" )
# BooleanOptionKey const correct_hbs_dihedrals ( "hbs_creator.correct_hbs_dihedrals" ) to be implemented if possible



if __name__ == '__main__':
	
	#option.add( hbs_creator.hbs_chain, "Chain from PDB to be mimicked. Default 'A'. Use letters." ).def("A")
	#option.add( hbs_creator.hbs_final_res, "Residue number of the final residue for mimicry. Default 1." ).def(1)
	#option.add( hbs_creator.hbs_length, "Number of residues to mimic. Default 12." ).def(12)
	#option.add( hbs_creator.final_repack, "Do a final repack. Default false" ).def(false)
	#option.add( hbs_creator.final_minimize, "Do a final minimization. Default false" ).def(false)
	#option.add( hbs_creator.final_mc, "Do a final monte carlo on hbs. Default false" ).def(false)
	#option.add( hbs_creator.correct_hbs_dihedrals, "Correct hbs dihedral to low energy well. Default false" ).def(false)
	
	# apply to pose.

def apply(pose):
	# create score function
	score_fxn = scoring.get_score_function()
	scoring.constraints.add_fa_constraints_from_cmdline_to_scorefxn(score_fxn)
	scoring.constraints.add_fa_constraints_from_cmdline_to_pose(pose)

	if score_fxn.has_zero_weight(atom_pair_constraint):
		score_fxn.set_weight(atom_pair_constraint, 0.1)
	if score_fxn.has_zero_weight( dihedral_constraint ):
		score_fxn.set_weight(dihedral_constraint, 0.1)
	if score_fxn.has_zero_weight( angle_constraint):
		score_fxn.set_weight(angle_constraint, 0.1)
	
	restype_set = chemical.ChemicalManager.get_instance().residue_type_set(core.chemical.FA_STANDARD)

	hbs_chain = option[hbs_creator.hbs_chain].value()[0]
	final_res = option[hbs_creator.hbs_final_res].value()
	pdb_info = pose.pdb_info()

	# shift the jump
	#core.kinematics.FoldTree f = pose.fold_tree()
	#f.slide_jump( 1, 1, pose.pdb_info().pdb2pose( hbs_chain, final_res+1 ) )
	#pose.fold_tree( f )

	Size patch_ros_num = 0

	for ii in xrange(pose.size()):
		i = ii + 1
		chn = pdb_info.chain(i)
		pdb_res_num = pdb_info.number(i)
		print "evaluating residue ", chn , " ", pdb_info.number(i)

		if chn != hbs_chain: continue
		# correct chain to be truncated and prepped

		# hbs pre is the smallest number of what we want to preserve
		while pdb_res_num < final_res:
			print "deleting residue ", pdb_res_num , " which was ", core.chemical.oneletter_code_from_aa(pose.aa(i))
			pose.delete_polymer_residue(i)
			pdb_res_num = pdb_info.number(i)

		if pdb_res_num > final_res + option[hbs_creator.hbs_length].value():
			#TR, "deleting residue ", pdb_res_num, std.endl
			while chn == hbs_chain && i <= pose.size():
				chn = pdb_info.chain(i)
				pose.delete_polymer_residue(i)

		if ( pdb_res_num == final_res ) patch_ros_num = i

	hbs_patcher = hbs.HbsPatcher(patch_ros_num)
	hbs_patcher.apply(pose)

	setup_pert_foldtree(pose)

	# presently the final residue in the pose is the terminal residue of the hbs
	# replace with terminal variant
	# AMW: This means that you have to clear other possible Cterm variant assignments
	# from this residue.
	conformation.Residue term(restype_set.get_residue_type_with_variant_added(
		restype_set.get_residue_type_with_variant_removed(
		pose.residue(pose.size()).type(),
		chemical.UPPER_TERMINUS_VARIANT),
		chemical.METHYLATED_CTERMINUS_VARIANT), True)

	term.set_all_chi(pose.residue(pose.size()).chi())
	#replace_res_post.mainchain_torsions(pose.residue(oop_post_pos_).mainchain_torsions())

	pose.replace_residue(pose.size(), term, True)
	conformation.idealize_position(pose.size(), pose.conformation())


	#pose.set_phi(i+3, -40)
	#pose.set_psi(i+3, -58)

	#pose.dump_pdb( "postdihedrals.pdb")

	pose.conformation().detect_bonds()
	pose.conformation().detect_pseudobonds()
	for i in xrange(pose.size()):
		pose.conformation().update_polymeric_connection(i + 1)

	pose.dump_pdb("postpseudobonds.pdb")

	# TINYMIN
	# create move map for minimization
	littlemm = kinematics.MoveMap()
	littlemm.set_bb(False)
	littlemm.set_chi(True)
	#mm.set_jump( 1, true )

	# create minimization mover
	littlemin = protocols.minimization_packing.MinMover(littlemm, score_fxn, option[ OptionKeys.run.min_type ].value(), 1, True)
	littlemin.apply(pose)

	if option[ hbs_creator.final_mc ].value():
		pert_sequence = moves.SequenceMover()
		pert_mc = moves.MonteCarlo(pose, score_fxn, 0.2)

		pert_pep_mm = kinematics.MoveMap()

		for ii in xrange(pose.size()):
			i = ii + 1
			if pdb_info.chain(i) == hbs_chain:
				if pose.residue_type(i).is_l_aa():
					print "setting small movable resid:", i
					#kdrew: commenting out because small mover fails randomly
					pert_pep_mm.set_bb(i)

		pert_pep_small = simple_moves.SmallMover(pert_pep_mm, 0.2, 1)
		pert_pep_small.angle_max('H', 2.0)
		pert_pep_small.angle_max('L', 2.0)
		pert_pep_small.angle_max('E', 2.0)

		pert_sequence.add_mover(pert_pep_small)

		#awatkins: add all hbs_pre positions to random small mover
		#TODO: I would PAY for understanding as to why this is so broken.
		#hbs.HbsRandomSmallMoverOP hpm( new hbs.HbsRandomSmallMover ( hbs_position, 2.0))#option[hbs_creator.hbs_length].value(), 2.0 ) )
		#moves.RepeatMoverOP pert_pep_repeat( new moves.RepeatMover( hpm, 1000 ) )
		#pert_sequence.add_mover( pert_pep_repeat )

		pert_trial = moves.TrialMover(pert_sequence, pert_mc)

		pert_trial.apply(pose)
		pert_mc.recover_low(pose)

	pose.dump_pdb("postmc.pdb")

	if option[hbs_creator.final_repack]:

		# create a task factory and task operations
		tf = TaskFactory()
		tf.push_back(operation.InitializeFromCommandline())

		#if ( ResourceManager.get_instance().has_option( packing.resfile ) ||  option[ packing.resfile ].user() ) {
		#	operation.ReadResfileOP rrop( new operation.ReadResfile() )
		#	rrop.default_filename()
		#	tf.push_back( rrop )
		#} else {
			#kdrew: do not do design, makes NATAA if res file is not specified
		rtrp = operation.RestrictToRepacking()
		tf.push_back(rtrp)
		#}


		# create a pack rotamers mover
		packer = protocols.minimization_packing.PackRotamersMover()
		packer.task_factory(tf)
		packer.score_function(score_fxn)

		packer.apply(pose)
	
	pose.dump_pdb("postrepack.pdb")

	if option[hbs_creator.final_minimize]:
		# create move map for minimization
		mm = kinematics.MoveMap()
		mm.set_chi(True)
		#mm.set_jump( 1, true )
		for ii in xrange(pose.size()):
			i = ii + 1
			if pose.pdb_info().chain(i) == hbs_chain:
				print "Can move ", i, " bb "
				mm.set_bb(i, True)
			else:
				mm.set_bb(i, False)
			
		print mm
		score_fxn.set_weight(atom_pair_constraint, 0.5)
		score_fxn.set_weight(dihedral_constraint, 0.5)
		score_fxn.set_weight(angle_constraint, 0.5)

		# create minimization mover
		minM = protocols.minimization_packing.MinMover(mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, True)
		#minM.cartesian( true )
		minM.apply(pose)

		score_fxn.set_weight(atom_pair_constraint, 1)
		score_fxn.set_weight(dihedral_constraint, 1)
		score_fxn.set_weight(angle_constraint, 1)
		minM.apply(pose)

# this only works for two chains and assumes the protein is first and the peptide is second
# inspired by protocols/docking/DockingProtocol.cc
def setup_pert_foldtree(pose):
	from kinematics import *

	# get current fold tree
	f = pose.fold_tree()
	f.clear()

	# get the start and end for both chains
	pro_start = pose.conformation().chain_begin(1)
	pro_end = pose.conformation().chain_end(1)
	pep_start = pose.conformation().chain_begin(2)
	pep_end = pose.conformation().chain_end(2)

	# get jump positions based on the center of mass of the chains
	dock_jump_pos_pro(core.pose.residue_center_of_mass(pose, pro_start, pro_end))
	dock_jump_pos_pep(core.pose.residue_center_of_mass(pose, pep_start, pep_end))

	# build fold tree
	jump_index(f.num_jump() + 1)
	#-1 is a magic number for PEPTIDE EDGE.  There is a constant defined with the fold tree that should have been used here.
	f.add_edge(pro_start, dock_jump_pos_pro, -1)
	f.add_edge(dock_jump_pos_pro, pro_end, -1)
	f.add_edge(pep_start, dock_jump_pos_pep, -1)
	f.add_edge(dock_jump_pos_pep, pep_end, -1)
	f.add_edge(dock_jump_pos_pro, dock_jump_pos_pep, jump_index)

	# set pose foldtree to foldtree we just created
	f.reorder(1)
	f.check_fold_tree()
	assert(f.check_fold_tree())

	print "AFTER:", f

	pose.fold_tree(f)
