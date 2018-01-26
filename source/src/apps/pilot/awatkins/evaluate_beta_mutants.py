#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   evaluate_beta_mutants.py
## @brief  aaa
## @author Andrew Watkins



def evaluate_interface(pose, score_fxn):
	v1 = scorefxn(pose)

	trans_mover = protocols.rigid.RigidBodyTransMover(pose, 1)
	trans_mover.step_size(1000)
	trans_mover.apply(pose)

	# repack min in unbound state.
	tf = TaskFactory()
	tf.push_back( core.pack.task.operation.InitializeFromCommandline())
	tf.push_back(operation.RestrictToRepacking())
	packer = protocols.minimization_packing.PackRotamersMover()
	packer.task_factory(tf)
	packer.score_function(score_fxn)
	packer.apply(pose)

	separate_min_mm = kinematics.MoveMap()
	separate_min_mm.set_bb(True)
	separate_min_mm.set_chi(True)
	separate_min_mm.set_jump(1, True)
	separate_min = minimization_packing.MinMover(separate_min_mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, True )
	separate_min.apply(pose)

	v2 = score_fxn(pose)
	return v1 - v2

def score_ensemble(pose, score_fxn, name, ntrials=3):
	frlx = protocols.relax.FastRelax(score_fxn)
	mm = kinematics.MoveMap()
	mm.set_bb(True)
	mm.set_chi(True)
	mm.set_jump(True)
	frlx.set_movemap(mm)

	# We actually want the multiple trials to occur on an unrelaxed pose
	# Otherwise we don't actually sample much.
	#frlx.apply( pose )

	score = evaluate_interface(pose, score_fxn)
	sum = 0
	print "Starting score of", name, "pose is:", score

	scores = []
	for ii in xrange(ntrials):
		copypose = Pose(pose)

		frlx.apply(copypose)

		scores.append(evaluate_interface(copypose, score_fxn))

		sum += scores[-1]
		if scores[-1] < score:
			score = scores[-1]
			copypose.dump_pdb(name)
		print "Trial", ii, "score", newscore

	print "Other stats: average was", sum/float(ntrials)

	# Av of best ten
	if ntrials > 10:
		print "Average of best 10: ", sum(sorted(scores)[:10])/10.0, std.endl
	return score

if __name__ == '__main__':

	# init command line options
	#you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
	#devel.init(argc, argv)

	
	# create score function
	score_fxn = core.scoring.get_score_function()
	score_fxn.set_weight( fa_rep, 0.3 )

	import sys
	Pose = core.import_pose.pose_from_file( pose, sys.argv[1], core.import_pose.PDB_file )
	bool is_1ycr = "1ycr" in sys.argv[1]
	bool is_mdm4 = "2z5t" in sys.argv[1]

	chm = rosetta.core.chemical.ChemicalManager.get_instance()
    restype_set = chm.residue_type_set('fa_standard').get_self_ptr()#.get_self_weak_ptr()

	# Hard coded because who cares.
	# For mdm2 based on 1ycr
	leu_pos = 96 if is_mdm4 else ( 88 if is_1ycr else 91 )
	trp_pos = leu_pos + 3
	phe_pos = leu_pos + 6

	leu = pose.residue(leu_pos)
	trp = pose.residue(trp_pos)
	phe = pose.residue(phe_pos)

	trp_mm = kinematics.MoveMap()
	trp_mm.set_bb(True, trp_pos)
	trp_mm.set_chi(True, trp_pos)
	trp_min = protocols.minimization_packing.MinMover(trp_mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.001, True)

	leu_mm = kinematics.MoveMap()
	leu_mm.set_bb(True, leu_pos)
	leu_mm.set_chi(True, leu_pos)
	leu_min = protocols.minimization_packing.MinMover(leu_mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.001, True)

	phe_mm = kinematics.MoveMap()
	phe_mm.set_bb(True, phe_pos)
	phe_mm.set_chi(True, phe_pos)
	phe_min = protocols.minimization_packing.MinMover(phe_mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.001, True)

	file = sys.argv[1]
	# This starting pose is the first
	print "Native."
	b53_8_score = score_ensemble(pose, score_fxn, "native_"+file)

	# Next, meta-trifluoromethyl.
	print "Mutant 1, meta-trifluoromethyl.", std.endl
	b53_12_pose = Pose(pose)
	tfm_rt = restype_set.name_map("B3F:B3F-CE1-trifluoromethylated")
	tfm = ResidueFactory.create_residue(tfm_rt, pose.residue(trp_pos), pose.conformation())
	core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(pose.residue(trp_pos), tfm, b53_12_pose.conformation())
	b53_12_pose.conformation().replace_residue(trp_pos, tfm, False)
	# Minimize the substituted residue.
	trp_min.apply(b53_12_pose)
	b53_12_score = score_ensemble(b53_12_pose, score_fxn, "meta_trifluoromethyl_phe_"+file)

	# Next, meta-trifluoromethyl.
	print "Mutant 1b, other meta-trifluoromethyl.", std.endl
	b53_12b_pose = Pose(pose)
	tfmb_rt = restype_set.name_map("B3F:B3F-CE2-trifluoromethylated")
	tfmb = ResidueFactory.create_residue(tfmb_rt, pose.residue(trp_pos), pose.conformation())
	core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(pose.residue(trp_pos), tfmb, b53_12b_pose.conformation())
	b53_12b_pose.conformation().replace_residue(trp_pos, tfmb, False)
	trp_min.apply(b53_12b_pose)
	b53_12b_score = score_ensemble(b53_12b_pose, score_fxn, "meta_trifluoromethyl_pheb_"+file)

	# Next, ch2-chloro trp
	print "Mutant 2, ch2-chloro trp."
	b53_13_pose = Pose(pose)
	ch2_rt = restype_set.name_map("B3W:B3W-CH2-chlorinated")
	ch2 = ResidueFactory.create_residue(ch2_rt, pose.residue(trp_pos), pose.conformation())
	core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(pose.residue(trp_pos), ch2, b53_13_pose.conformation())
	b53_13_pose.conformation().replace_residue(trp_pos, ch2, False)
	trp_min.apply(b53_13_pose)
	b53_13_score = score_ensemble(b53_13_pose, score_fxn, "chloro_trp_"+file)

	# Next, phe
	print "Mutant 3, phe."
	b53_14_pose = Pose(pose)
	phe2_rt = restype_set.name_map("B3F")
	phe2 = ResidueFactory.create_residue(phe2_rt, pose.residue(trp_pos), pose.conformation())
	core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(pose.residue(trp_pos), phe2, b53_14_pose.conformation())
	b53_14_pose.conformation().replace_residue(trp_pos, phe2, False)
	trp_min.apply(b53_14_pose)
	b53_14_score = score_ensemble(b53_14_pose, score_fxn, "phe_"+file)

	# Next, meta-chloro
	print "Mutant 4, m-chlorophe."
	b53_15_pose = Pose(pose)
	mcl_rt = restype_set.name_map("B3F:B3F-CE1-chlorinated")
	mcl = ResidueFactory.create_residue(mcl_rt, pose.residue(trp_pos), pose.conformation())
	core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(pose.residue(trp_pos), mcl, b53_15_pose.conformation())
	b53_15_pose.conformation().replace_residue(trp_pos, mcl, False)
	trp_min.apply(b53_15_pose)
	b53_15_score = score_ensemble(b53_15_pose, score_fxn, "meta_chloro_phe_"+file)

	# Next, meta para dichloro
	print "Mutant 5, m,p dichloro phe."
	b53_16_pose = Pose(pose)
	mpcl_rt = restype_set.name_map("B3F:B3F-CE1-chlorinated:B3F-CZ-chlorinated")
	mpcl = ResidueFactory.create_residue(mpcl_rt, pose.residue(trp_pos), pose.conformation())
	core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(pose.residue(trp_pos), mpcl, b53_16_pose.conformation())
	b53_16_pose.conformation().replace_residue(trp_pos, mpcl, False)
	trp_min.apply(b53_16_pose)
	b53_16_score = score_ensemble(b53_16_pose, score_fxn, "meta_para_dichloro_phe_"+file)

	# Next, leu . ile from 12
	print "Mutant 6, ile and also m-trifluoromethyl."
	b53_17_pose = Pose(b53_12_pose)
	ile_rt = restype_set.name_map( "B3I" )
	ile = ResidueFactory.create_residue( ile_rt, pose.residue(leu_pos), pose.conformation() )
	core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(pose.residue(leu_pos), ile, b53_17_pose.conformation())
	b53_17_pose.conformation().replace_residue(leu_pos, ile, False)
	leu_min.apply(b53_17_pose)
	b53_17_score = score_ensemble(b53_17_pose, score_fxn, "ile_"+file)
	
	# Next, phe . 4-fl phe from 16
	print "Mutant 7, m,p dichloro phe and also 4-fl phe."
	b53_18_pose = Pose(b53_16_pose)
	fp_rt = restype_set.name_map("B3F:B3F-CZ-fluorinated")
	fp = ResidueFactory.create_residue(fp_rt, pose.residue(phe_pos), pose.conformation())
	core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(pose.residue(phe_pos), fp, b53_18_pose.conformation())
	b53_18_pose.conformation().replace_residue(phe_pos, fp, False)
	phe_min.apply(b53_18_pose)
	b53_18_score = score_ensemble(b53_18_pose, score_fxn, "para_fluoro_phe_"+file)
	
	print "Starting score beta-53-8 for mdm2: ", b53_8_score
	print "meta trifluoromethyl phe beta-53-12 for mdm2: ", b53_12_scor
	print "meta trifluoromethyl phe beta-53-12 for mdm2: ", b53_12b_score
	print "CH2-chloro trp beta-53-13 for mdm2: ", b53_13_score
	print "Phe beta-53-14 for mdm2: ", b53_14_score
	print "meta chloro phe beta-53-15 for mdm2: ", b53_15_score
	print "meta para dichloro phe beta-53-16 for mdm2: ", b53_16_score
	print "meta trifluoromethyl phe ile beta-53-17 for mdm2: ", b53_17_score
	print "para fluoro phe meta para dichloro phe beta-53-18 for mdm2: ", b53_18_score
