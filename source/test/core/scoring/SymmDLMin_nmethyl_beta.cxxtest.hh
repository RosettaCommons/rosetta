// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file test/core/scoring/SymmDLMin_nmethyl_beta.cxxtest.hh
/// @brief Unit tests for minimiziation of mirror-image poses with the --symmetric_gly_tables option and the beta scorefunction,
/// with some N-methylated residues thrown in.
/// @detials Mirror image poses (with D- and L-amino acids swapped) should minimize identically with this option.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pdb1rpb.hh>
#include <test/core/init_util.hh>
#include <test/util/symmetry_funcs.hh>

// Unit headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/chemical/ResidueType.hh>
#include <protocols/simple_moves/MutateResidue.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>

//Minimizer
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/MoveMap.hh>

#include <basic/Tracer.hh> // AUTO IWYU For Tracer

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;


static basic::Tracer TR("core.scoring.SymmDLMin_nmethyl_beta_Tests.cxxtest");

class SymmDLMin_nmethyl_beta_Tests : public CxxTest::TestSuite {

public:

	void setUp() {
		//core_init();
		core_init_with_additional_options( "-beta_nov16 -score:weights beta_nov16.wts -symmetric_gly_tables true -write_all_connect_info -connect_info_cutoff 0.0" );
	}

	void tearDown() {
	}

	/// @brief Run the minimizer on the pose.
	///
	void do_minimization( core::pose::Pose &pose, core::scoring::ScoreFunctionOP sfxn, bool const cartesian ) {
		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb(true);
		mm->set_chi(true);
		if ( cartesian ) {
			core::optimization::CartesianMinimizer minimizer;
			core::optimization::MinimizerOptionsOP min_options( new core::optimization::MinimizerOptions( "linmin", 10.0, true, false, false ) );
			minimizer.run( pose, *mm, *sfxn, *min_options );
		} else {
			core::optimization::AtomTreeMinimizer minimizer;
			core::optimization::MinimizerOptionsOP min_options( new core::optimization::MinimizerOptions( "linmin", 10.0, true, false, false ) );
			minimizer.run( pose, *mm, *sfxn, *min_options );
		}
	}

	/// @brief Are two angles within a threshhold of one another?
	///
	bool within_thresh( core::Real const &val1, core::Real const &val2, core::Real const &thresh ) {
		if ( std::abs( val1 - val2 ) <= thresh ) return true;
		core::Real val1prime = val1;
		core::Real val2prime = val2;
		if ( val1<val2 ) {
			val1prime += 360.0;
		} else {
			val2prime += 360.0;
		}
		return ( std::abs( val1prime - val2prime ) <= thresh );
	}

	/// @brief Construct a few L-amino acid poses, and confirm that mirror-image conformations
	/// score identically with a given scorefunction.
	void mirror_pose_test( core::scoring::ScoreFunctionOP sfxn, bool const cartesian ) {
		core::pose::Pose pose = pdb1rpb_pose();

		// Add peptoids:
		for ( core::Size i=2, imax=pose.total_residue(); i<imax; i+=3 ) { //Make every first and second residue (of every three) a peptoid
			protocols::simple_moves::MutateResidue mutres( i, "ALA:N_Methylation" );
			protocols::simple_moves::MutateResidue mutres2( i-1, "TRP:N_Methylation" );
			mutres.set_update_polymer_dependent( true );
			mutres2.set_update_polymer_dependent( true );
			if ( !pose.residue_type(i).is_disulfide_bonded() ) {
				mutres.apply(pose);
				//pose.set_chi( 1, i, -93.73 );
			}

			if ( !pose.residue_type(i-1).is_disulfide_bonded() ) {
				mutres2.apply(pose);
				pose.set_chi( 1, i-1, -56.9 );
				pose.set_chi( 2, i-1, -88.1 );
			}
		}
		pose.update_residue_neighbors();

		core::pose::Pose pose2 = pose;

		flip_residues(pose2);
		mirror_pose(pose, pose2);
		pose2.update_residue_neighbors();

		(*sfxn)(pose);
		(*sfxn)(pose2);
		//for(core::Size j=1; j<=100; ++j) {
		do_minimization(pose, sfxn, cartesian);
		do_minimization(pose2, sfxn, cartesian);
		//}

		(*sfxn)(pose);
		(*sfxn)(pose2);

		TS_ASSERT_DELTA(pose.energies().total_energy(), pose2.energies().total_energy(), std::abs( std::max(pose.energies().total_energy(), pose2.energies().total_energy())/100.0 ) );
		for ( core::Size ir=1, irmax=pose.size(); ir<=irmax; ++ir ) {
			TR << ir << "\tphiL=" << pose.phi(ir)   << "\tphiD=" << pose2.phi(ir)   << std::endl;
			TR << ir << "\tpsiL=" << pose.psi(ir)   << "\tpsiD=" << pose2.psi(ir)   << std::endl;
			TR << ir << "\tomgL=" << pose.omega(ir) << "\tomgD=" << pose2.omega(ir) << std::endl;
			TS_ASSERT( within_thresh( pose.phi(ir), -1.0*pose2.phi(ir), 0.015 ) );
			TS_ASSERT( within_thresh( pose.psi(ir), -1.0*pose2.psi(ir), 0.015 ) );
			TS_ASSERT( within_thresh( pose.omega(ir), -1.0*pose2.omega(ir), 0.015 ) );
			for ( core::Size ichi=1, ichimax=pose.residue(ir).nchi(); ichi<=ichimax; ++ichi ) {
				TR << ir << "\tchi" << ichi << "L=" << pose.chi(ichi,ir) << "\tchi" << ichi << "D=" << pose2.chi(ichi,ir) << std::endl;
				TS_ASSERT( within_thresh( pose.chi(ichi, ir), -1.0*pose2.chi(ichi, ir), 0.015 ) );
			}
		}
		TR.flush();
		//pose.dump_scored_pdb( "Ltemp.pdb", *sfxn );
		//core::pose::Pose pose3(pose2);
		//mirror_pose(pose2,pose3);
		//pose3.dump_scored_pdb( "Dtemp.pdb", *sfxn );
	}

	/// @brief Construct a few L-amino acid poses, and confirm that mirror-image conformations
	/// score identically with a given scorefunction.
	/// @details This version ensures that there's a GLY-PRO sequence.
	void mirror_pose_test2( core::scoring::ScoreFunctionOP sfxn, bool const cartesian ) {
		core::pose::Pose pose = pdb1rpb_pose();

		//Mutate residue 4 to 601.  Since residue 3 is a gly, this makes a nice gly-601.
		protocols::simple_moves::MutateResidue mutres;
		mutres.set_res_name("ALA:N_Methylation");
		mutres.set_target(4);
		mutres.apply(pose);

		//Mutate residues 2 and 3 to peptoids.
		protocols::simple_moves::MutateResidue mutres2;
		protocols::simple_moves::MutateResidue mutres3;
		mutres2.set_res_name("TRP:N_Methylation");
		mutres3.set_res_name("ALA:N_Methylation");
		mutres2.set_update_polymer_dependent(true);
		mutres3.set_update_polymer_dependent(true);
		mutres2.set_target(2);
		mutres3.set_target(3);
		mutres2.apply(pose);
		mutres3.apply(pose);

		//Set the sidechains to something
		pose.set_chi( 1, 2, 43.2 );
		pose.set_chi( 2, 2, -82.2 );
		pose.update_residue_neighbors();

		core::pose::Pose pose2 = pose;

		flip_residues(pose2);
		mirror_pose(pose, pose2);
		pose2.update_residue_neighbors();

		(*sfxn)(pose);
		(*sfxn)(pose2);
		//for(core::Size j=1; j<=100; ++j) {
		do_minimization(pose, sfxn, cartesian);
		do_minimization(pose2, sfxn, cartesian);
		//}

		(*sfxn)(pose);
		(*sfxn)(pose2);

		TS_ASSERT_DELTA(pose.energies().total_energy(), pose2.energies().total_energy(), std::abs( std::max(pose.energies().total_energy(), pose2.energies().total_energy())/100.0 ) );
		for ( core::Size ir=1, irmax=pose.size(); ir<=irmax; ++ir ) {
			TR << ir << "\tphiL=" << pose.phi(ir)   << "\tphiD=" << pose2.phi(ir)   << std::endl;
			TR << ir << "\tpsiL=" << pose.psi(ir)   << "\tpsiD=" << pose2.psi(ir)   << std::endl;
			TR << ir << "\tomgL=" << pose.omega(ir) << "\tomgD=" << pose2.omega(ir) << std::endl;
			TS_ASSERT( within_thresh( pose.phi(ir), -1.0*pose2.phi(ir), 0.015 ) );
			TS_ASSERT( within_thresh( pose.psi(ir), -1.0*pose2.psi(ir), 0.015 ) );
			TS_ASSERT( within_thresh( pose.omega(ir), -1.0*pose2.omega(ir), 0.015 ) );
			for ( core::Size ichi=1, ichimax=pose.residue(ir).nchi(); ichi<=ichimax; ++ichi ) {
				TR << ir << "\tchi" << ichi << "L=" << pose.chi(ichi,ir) << "\tchi" << ichi << "D=" << pose2.chi(ichi,ir) << std::endl;
				TS_ASSERT( within_thresh( pose.chi(ichi, ir), -1.0*pose2.chi(ichi, ir), 0.015 ) );
			}
		}
		TR.flush();
		//pose.dump_scored_pdb( "Ltemp.pdb", *sfxn );
		//core::pose::Pose pose3(pose2);
		//mirror_pose(pose2,pose3);
		//pose3.dump_scored_pdb( "Dtemp.pdb", *sfxn );
	}

	/// @brief Tests symmetric scoring with the cart_bonded scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_cart_bonded() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::cart_bonded, 1.0 );
		TR << "Testing cart_bonded score term." << std::endl;
		mirror_pose_test(scorefxn, true);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_atr scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_fa_atr() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_atr, 1.0 );
		TR << "Testing fa_atr score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_fa_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_rep, 1.0 );
		TR << "Testing fa_rep score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_intra_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_fa_intra_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_rep, 1.0 );
		TR << "Testing fa_intra_rep score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_sol scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_fa_sol() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_sol score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_intra_sol_xover4 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_fa_intra_sol_xover4() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_sol_xover4, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_intra_sol_xover4 score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the lk_ball_wtd scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_lk_ball_wtd() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::lk_ball_wtd, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing lk_ball_wtd score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_elec scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_fa_elec() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_elec, 1.0 );
		TR << "Testing fa_elec score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the pro_close scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_pro_close() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::pro_close, 1.0 );
		TR << "Testing pro_close score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the dslf_fa13 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_dslf_fa13() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::dslf_fa13, 1.0 );
		TR << "Testing dslf_fa13 score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the hbonds scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_hbonds() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::hbond_sr_bb, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_lr_bb, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_sc, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_bb_sc, 1.0 );
		TR << "Testing hbonds score terms." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_dun scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_fa_dun() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_dun, 1.0 );
		TR << "Testing fa_dun score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the omega scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_omega() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::omega, 1.0 );
		TR << "Testing omega score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the rama scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_rama() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama, 1.0 );
		TR << "Testing rama score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the rama_prepro scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_rama_prepro() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama_prepro, 1.0 );
		TR << "Testing rama_prepro score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the rama_prepro scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_rama_prepro2() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama_prepro, 1.0 );
		TR << "Testing rama_prepro score term with a gly-pro sequence." << std::endl;
		mirror_pose_test2(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the p_aa_pp scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_p_aa_pp() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::p_aa_pp, 1.0 );
		TR << "Testing p_aa_pp score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the yhh_planarity scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_yhh_planarity() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::yhh_planarity, 1.0 );
		TR << "Testing yhh_planarity score term." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the full beta_nov15 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_beta_sfxn() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->add_weights_from_file("beta_nov16.wts");
		TR << "Testing full beta_nov16 score function." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring with the full scorefunction that's the current default (whatever that is).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_min_current_default_scorefxn() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
		TR << "Testing full default score function." << std::endl;
		mirror_pose_test(scorefxn, false);
		return;
	}

};
