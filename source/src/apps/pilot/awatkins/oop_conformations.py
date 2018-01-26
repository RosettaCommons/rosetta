#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   best_beta_backbones.py
## @brief  Creates an OOP dimer and toys around with its dihedrals.
## @author Andrew Watkins

#namespace oop_conformations {
#// pert options
#RealOptionKey const min_threshold ( "oop_conformations.min_threshold" );
#RealOptionKey const dihedral_start ( "oop_conformations.dihedral_start" );
#RealOptionKey const dihedral_end ( "oop_conformations.dihedral_end" );
#BooleanOptionKey const oop_optimize( "oop_conformations.oop_optimize" );
#IntegerOptionKey const length( "oop_conformations.length" );#

def idealize(pose, scorefxn):

	pert_sequence = moves.SequenceMover()
	pert_mc = moves.MonteCarlo(pose, scorefxn, 0.2)
	pert_pep_mm = kinematics.MoveMap()

	oop_pre_positions = []
	for im1 in xrange(pose.size()):
		ii = im1 + 1
		if pose.residue(ii).has_variant_type(chemical.OOP_PRE):
			oop_pre_positions.append(ii)
		else:
			pert_pep_mm.set_bb(ii, True )

	pert_pep_small = simple_moves.SmallMover(pert_pep_mm, 0.2, 1 )
	pert_pep_small.angle_max('H', 2.0)
	pert_pep_small.angle_max('L', 2.0)
	pert_pep_small.angle_max('E', 2.0)

	pert_sequence.add_mover( pert_pep_small );

	if len(oop_pre_positions) > 0:
		opm = simple_moves.oop.OopRandomSmallMover(oop_pre_positions, 2.0)
		pert_pep_repeat = moves.RepeatMover(opm, len(oop_pre_positions) * 50)
		pert_sequence.add_mover(pert_pep_repeat)
	pert_trial = moves.TrialMover(pert_sequence, pert_mc)

	pert_trial.apply(pose)
	pert_mc.recover_low(pose)

if __name__ == '__main__':

	#option.add( oop_conformations.min_threshold, "Minimization threshold" ).def(100);
	#option.add( oop_conformations.dihedral_start, "Start dihedral" ).def(30);
	#option.add( oop_conformations.dihedral_end, "End dihedral." ).def(360);
	#option.add( oop_conformations.oop_optimize, "Do a final conformational sampling. Default true" ).def(true);
	#option.add( oop_conformations.length, "Number residues" ).def(4);

	#	devel.init( argc, argv );

	Pose ala_pose;

	# Make an oop.
	length = option[oop_conformations.length].value();
	ala_pose = pyrosetta.pose_from_sequence("A"*length)

	plus_pos = core.select.residue_selector.ResidueIndexSelector()
	for i in xrange(length): plus_pos.append_index(i + 1)

	OC_mover = protocols.ncbb.oop.OopCreatorMover(plus_pos, 0, 0, 0, 0, 0, 0, False, True, False, True );
	
	mvmp = kinematics.MoveMap()
	scorefxn_no_hbond = ScoreFunctionFactory.create_score_function("mm_std_no_hbond")
	scorefxn = ScoreFunctionFactory.create_score_function("mm_std")

	if scorefxn.has_zero_weight(core.scoring.atom_pair_constraint):
		scorefxn.set_weight(core.scoring.atom_pair_constraint, 1.0)

	mvmp.set_chi(1, True)
	mvmp.set_bb( True );
	mnmvr = protocols.minimization_packing.MinMover(mvmp, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, True)
	mnmvr_no_hbond = protocols.minimization_packing.MinMover(mvmp, scorefxn_no_hbond, "lbfgs_armijo_nonmonotone", 0.0001, True)
	#mnmvr.cartesian( true );
	#mnmvr_no_hbond.cartesian( true );

	# MUST minimize initially or pose is ridiculous
	#mnmvr_no_hbond.apply( pose );
	best_angles = [0 for _ in xrange(8)]

	Real lowest_energy = 10000;
	Pose best_pose = ala_pose;

	# Do 10,000 samples, arbitrarily!
	for nn in xrange(10000):
		Pose pose( ala_pose );

		angles = [360 for _ in xrange(8)]
		for resi in xrange(pose.size()):
			while angles[resi * 2 ] >= 180 or angles[resi*2] <= -180:
				angles[resi * 2 ] = numeric.random.rg().uniform() * 359.8 - 179.9;
			while angles[resi*2 + 1] >= 180 or angles[resi*2 + 1] <= -180:
				angles[resi * 2 + 1 ] = numeric.random.rg().uniform() * 359.8 - 179.9;
			pose.set_phi( resi, angles[resi * 2 - 1] );
			pose.set_psi( resi, angles[resi * 2 ] );
		OC_mover.apply( pose );

		print "\nTrial", n
		print "Best energy so far is ", lowest_energy, "."

		# score the pose
		orig_ener = scorefxn(pose);
		print angles
		for im1 in xrange(pose.size()):
			print "({phi}, {psi}), ".format(phi=pose.phi(im1+1), psi=pose.psi(im1+1)),
		print " and energy", orig_ener,

		mnmvr_no_hbond.apply( pose );
		print "minimizes to", scorefxn(pose)

		#if ( option[oop_conformations.oop_optimize].value() ) {
		idealize( pose, scorefxn );
		print "OOP moves to ", scorefxn(pose)
		
		if scorefxn(pose) < 10:
			mnmvr.apply( pose );
			print "More minimization to ", scorefxn(pose)
		


		energy = scorefxn(pose)
		if energy < lowest_energy:
			best_pose = pose;
			lowest_energy = energy;
			real_phi = [pose.phi(resi+1) for resi in xrange(pose.size())]
			real_psi = [pose.psi(resi+1) for resi in xrange(pose.size())]

			print "new min ", lowest_energy," found at ",
			for phi, psi in zip(real_phi, real_psi):
				print "({phi}, {psi})".format(phi=phi, psi=psi),
			print ""
			for resi in xrange(pose.size()):
				best_angles[ 2 * resi     ] = real_phi[ resi ]
				best_angles[ 2 * resi + 1 ] = real_psi[ resi ]

	print "Final lowest energy is", lowest_energy. "found at", best_angles

	best_pose.dump_pdb("final_no_final_opt.pdb")
	idealize(best_pose, scorefxn)
	mnmvr.apply(best_pose )
	best_pose.dump_pdb("final_final_opt.pdb" )
