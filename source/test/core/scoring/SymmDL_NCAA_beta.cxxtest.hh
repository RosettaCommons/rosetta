// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file test/core/scoring/SymmDL_NCAA_beta.cxxtest.hh
/// @brief Unit tests for symmetric pose scoring with the --symmetric_gly_tables option and noncanonicals, with the beta scorefunction.
/// @detials Mirror image poses (with D- and L-amino acids swapped) should score identically with this option.
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
#include <protocols/simple_moves/MutateResidue.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>

//Minimizer

#include <basic/Tracer.hh> // AUTO IWYU For Tracer

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;


static basic::Tracer TR("core.scoring.SymmDL_NCAA_betaTests.cxxtest");

class SymmDL_NCAA_betaTests : public CxxTest::TestSuite {

public:

	void setUp() {
		//core_init();
		core_init_with_additional_options( "-beta_nov16 -symmetric_gly_tables true -write_all_connect_info -connect_info_cutoff 0.0" );
	}

	void tearDown() {
	}

	/// @brief Construct a few L-amino acid poses, and confirm that mirror-image conformations
	/// score identically with a given scorefunction.
	void mirror_pose_test( core::scoring::ScoreFunctionOP sfxn ) {
		core::pose::Pose pose = pdb1rpb_pose();

		// Mutate a few residues to noncanonicals
		protocols::simple_moves::MutateResidue mutres1;
		mutres1.set_res_name("ORN");
		mutres1.set_target(4);
		mutres1.apply(pose);
		protocols::simple_moves::MutateResidue mutres2;
		mutres2.set_res_name("DAB");
		mutres2.set_target(5);
		mutres2.apply(pose);
		protocols::simple_moves::MutateResidue mutres3;
		mutres3.set_res_name("DPP");
		mutres3.set_target(2);
		mutres3.apply(pose);

		// Set sidechains:
		pose.set_chi( 1, 4, 145.3 );
		pose.set_chi( 2, 4, -163.8 );
		pose.set_chi( 3, 4, 66.7 );
		pose.set_chi( 1, 5, -67.3 );
		pose.set_chi( 2, 5, 170.3 );
		pose.set_chi( 2, 5, 160.5 );
		pose.update_residue_neighbors();

		core::pose::Pose pose2 = pose;
		flip_residues(pose2);
		mirror_pose(pose, pose2);
		pose2.update_residue_neighbors();

		(*sfxn)(pose);
		(*sfxn)(pose2);
		TS_ASSERT_DELTA(pose.energies().total_energy(), pose2.energies().total_energy(), std::abs( std::max(pose.energies().total_energy(), pose2.energies().total_energy())/1000.0 ) );
		//pose.dump_scored_pdb( "Ltemp.pdb", *sfxn );
		//pose2.dump_scored_pdb( "Dtemp.pdb", *sfxn );
	}


	/// @brief Construct a few L-amino acid poses, and confirm that mirror-image conformations
	/// score identically with a given scorefunction.
	/// @details This version ensures that there's a glycine-proline pair.
	void mirror_pose_test2( core::scoring::ScoreFunctionOP sfxn ) {
		core::pose::Pose pose = pdb1rpb_pose();

		//Mutate residue 4 to proline and residue 3 to DORN.  This makes a nice NCAA-pro test.
		protocols::simple_moves::MutateResidue mutres;
		mutres.set_res_name("PRO");
		mutres.set_target(4);
		mutres.apply(pose);
		protocols::simple_moves::MutateResidue mutres2;
		mutres2.set_res_name("DORN");
		mutres2.set_target(3);
		mutres2.apply(pose);

		core::pose::Pose pose2 = pose;

		flip_residues(pose2);
		mirror_pose(pose, pose2);
		pose2.update_residue_neighbors();

		(*sfxn)(pose);
		(*sfxn)(pose2);
		TS_ASSERT_DELTA(pose.energies().total_energy(), pose2.energies().total_energy(), std::abs( std::max(pose.energies().total_energy(), pose2.energies().total_energy())/1000.0 ) );
		//pose.dump_scored_pdb( "Ltemp.pdb", *sfxn );
		//pose2.dump_scored_pdb( "Dtemp.pdb", *sfxn );
	}

	/// @brief Tests symmetric scoring with the cart_bonded scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_cart_bonded() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::cart_bonded, 1.0 );
		TR << "Testing cart_bonded score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_atr scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_fa_atr() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::fa_atr, 1.0 );
		TR << "Testing fa_atr score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_fa_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::fa_rep, 1.0 );
		TR << "Testing fa_rep score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_intra_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_fa_intra_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::fa_intra_rep, 1.0 );
		TR << "Testing fa_intra_rep score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_sol scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_fa_sol() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_sol score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_elec scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_fa_elec() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::fa_elec, 1.0 );
		TR << "Testing fa_elec score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the pro_close scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_pro_close() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::pro_close, 1.0 );
		TR << "Testing pro_close score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the dslf_fa13 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_dslf_fa13() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::dslf_fa13, 1.0 );
		TR << "Testing dslf_fa13 score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the hbonds scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_hbonds() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::hbond_sr_bb, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_lr_bb, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_sc, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_bb_sc, 1.0 );
		TR << "Testing hbonds score terms." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_dun scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_fa_dun() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::fa_dun, 1.0 );
		TR << "Testing fa_dun score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the omega scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_omega() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::omega, 1.0 );
		TR << "Testing omega score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the rama_prepro scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_rama_prepro() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::rama_prepro, 1.0 );
		TR << "Testing rama_prepro score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the rama_prepro scorefunction and a gly-pro pair.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_rama_prepro2() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::rama_prepro, 1.0 );
		TR << "Testing rama_prepro score term with a gly-pro pair." << std::endl;
		mirror_pose_test2(scorefxn);
		return;
	}


	/// @brief Tests symmetric scoring with the p_aa_pp scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_p_aa_pp() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::p_aa_pp, 1.0 );
		TR << "Testing p_aa_pp score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the yhh_planarity scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_yhh_planarity() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->set_weight( core::scoring::yhh_planarity, 1.0 );
		TR << "Testing yhh_planarity score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the full talaris2014 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_beta() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		scorefxn->add_weights_from_file("beta_nov16.wts");
		TR << "Testing full beta_nov16 score function." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the full default scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_default_beta_scorefxn() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
		TR << "Testing full beta_nov16 score function." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

};
