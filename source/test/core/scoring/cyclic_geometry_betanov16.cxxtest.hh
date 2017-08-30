// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file test/core/scoring/cyclic_geometry.cxxtest.hh
/// @brief Unit tests for cyclic peptide pose scoring.
/// @detials Cyclic permutations should score identically.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>
#include <test/util/symmetry_funcs.hh>

// Unit headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/util.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/chemical/AA.hh>

#include <test/core/scoring/cyclic_geometry_headers.hh>

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;


static THREAD_LOCAL basic::Tracer TR("core.scoring.CyclicGeometry_betanov16_Tests.cxxtest");

class CyclicGeometry_betanov16_Tests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options( "-beta_nov16 -score:weights beta_nov16.wts -symmetric_gly_tables true -write_all_connect_info -connect_info_cutoff 0.0" );

		// Pull in the cyclic peptide pose (9 residues):
		core::pose::PoseOP initial_pose( new core::pose::Pose );
		core::import_pose::pose_from_file( *initial_pose, "core/scoring/cyclic_peptide.pdb" , core::import_pose::PDB_file );

		// Strip off termini and connect the ends:
		core::pose::remove_variant_type_from_pose_residue( *initial_pose, core::chemical::LOWER_TERMINUS_VARIANT, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose, core::chemical::CUTPOINT_LOWER, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose, core::chemical::CUTPOINT_UPPER, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose, core::chemical::UPPER_TERMINUS_VARIANT, 9 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose, core::chemical::CUTPOINT_LOWER, 9 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose, core::chemical::CUTPOINT_UPPER, 9 );
		initial_pose->conformation().declare_chemical_bond(1, "N", 9, "C");
		for ( core::Size i(1); i<=9; ++i ) {
			initial_pose->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(i);
		}

		poses_.clear();
		mirror_poses_.clear();

		poses_.push_back(initial_pose);
		mirror_poses_.push_back( mirror_pose_with_disulfides( poses_[1] ) );
		for ( core::Size i=1; i<=8; ++i ) {
			poses_.push_back( permute( poses_[i] ) );
			mirror_poses_.push_back( mirror_pose_with_disulfides( poses_[i+1] ) );
		}
	}

	void tearDown() {
	}

	/// @brief Tests cyclic permutation scoring with the cart_bonded scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_cart_bonded() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::cart_bonded, 1.0 );
		TR << "Testing cart_bonded score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_, false);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_atr scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_atr() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_atr, 1.0 );
		TR << "Testing fa_atr score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_rep, 1.0 );
		TR << "Testing fa_rep score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_intra_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_intra_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_rep, 1.0 );
		TR << "Testing fa_intra_rep score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_sol scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_sol() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_sol score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_intra_atr_xover4 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_intra_atr_xover4() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_atr_xover4, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_intra_atr_xover4 score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_intra_rep_xover4 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_intra_rep_xover4() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_rep_xover4, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_intra_rep_xover4 score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_intra_sol_xover4 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_intra_sol_xover4() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_sol_xover4, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_intra_sol_xover4 score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the lk_ball scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_lk_ball() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::lk_ball, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing lk_ball score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the lk_ball_iso scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_lk_ball_iso() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::lk_ball_iso, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing lk_ball_iso score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the lk_ball_bridge scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_lk_ball_bridge() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::lk_ball_bridge, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing lk_ball_bridge score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the lk_ball_bridge_uncpl scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_lk_ball_bridge_uncpl() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::lk_ball_bridge_uncpl, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing lk_ball_bridge_uncpl score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_elec scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_elec() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_elec, 1.0 );
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_intra_elec scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_intra_elec() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_elec, 1.0 );
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the pro_close scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_pro_close() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::pro_close, 1.0 );
		TR << "Testing pro_close score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the dslf_fa13 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_dslf_fa13() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::dslf_fa13, 1.0 );
		TR << "Testing dslf_fa13 score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the hbonds scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_hbonds() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::hbond_sr_bb, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_lr_bb, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_sc, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_bb_sc, 1.0 );
		TR << "Testing hbonds score terms." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_dun_rot scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_dun_rot() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_dun_rot, 1.0 );
		TR << "Testing fa_dun_rot score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_dun_dev scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_dun_dev() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_dun_dev, 1.0 );
		TR << "Testing fa_dun_dev score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_dun_semi scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_dun_semi() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_dun_semi, 1.0 );
		TR << "Testing fa_dun_semi score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the hxl_tors scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_hxl_tors() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::hxl_tors, 1.0 );
		TR << "Testing hxl_tors score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the omega scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_omega() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::omega, 1.0 );
		TR << "Testing omega score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the rama scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_rama() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama, 1.0 );
		TR << "Testing rama score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the rama_prepro scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_rama_prepro() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama_prepro, 1.0 );
		TR << "Testing rama_prepro score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the p_aa_pp scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_p_aa_pp() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::p_aa_pp, 1.0 );
		TR << "Testing p_aa_pp score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the full beta_nov16 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_beta_nov16() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->add_weights_from_file("beta_nov16.wts");
		TR << "Testing full beta_nov16 score function." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the full scorefunction that's the current default (whatever that is).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_current_default_scorefxn() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
		TR << "Testing full default score function." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_, mirror_poses_);
		return;
	}

private:
	utility::vector1 < core::pose::PoseOP > poses_;
	utility::vector1 < core::pose::PoseOP > mirror_poses_;


};
