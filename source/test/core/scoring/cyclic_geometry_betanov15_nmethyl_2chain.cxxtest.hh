// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file test/core/scoring/cyclic_geometry_betanov15_nmethyl_2chain.cxxtest.hh
/// @brief Unit tests for cyclic peptide pose scoring with 2 chains in the peptide and the beta_nov15 score function, with N-methylation.
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
#include <core/pose/variant_util.hh>
#include <core/chemical/AA.hh>

//Protocols headers
#include <protocols/simple_moves/MutateResidue.hh>

#include <test/core/scoring/cyclic_geometry_headers.hh>

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;


static THREAD_LOCAL basic::Tracer TR("core.scoring.CyclicGeometry_betanov15_Tests.cxxtest");

class CyclicGeometry_nmethyl_betanov15_TwoChainTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options( "-symmetric_gly_tables true -write_all_connect_info -connect_info_cutoff 0.0 -score:weights beta_nov15.wts" );

		// Pull in the two-chain cyclic peptide pose (12 + 12 = 24 residues)
		core::pose::PoseOP initial_pose_2chain( new core::pose::Pose );
		core::import_pose::pose_from_file( *initial_pose_2chain, "core/scoring/twochain_cycpep.pdb", core::import_pose::PDB_file );

		// Strip off termini and connect the ends:
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::LOWER_TERMINUS_VARIANT, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_LOWER, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_UPPER, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::UPPER_TERMINUS_VARIANT, 12 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_LOWER, 12 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_UPPER, 12 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::LOWER_TERMINUS_VARIANT, 13 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_LOWER, 13 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_UPPER, 13 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::UPPER_TERMINUS_VARIANT, 24 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_LOWER, 24 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_UPPER, 24 );
		initial_pose_2chain->conformation().declare_chemical_bond(1, "N", 24, "C");
		initial_pose_2chain->conformation().declare_chemical_bond(13, "N", 12, "C");

		remove_disulfides(initial_pose_2chain);
		form_disulfides(initial_pose_2chain);
		for ( core::Size ir=1, irmax=initial_pose_2chain->size(); ir<=irmax; ++ir ) {
			initial_pose_2chain->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(ir);
		}

		// Add N-methylation:
		protocols::simple_moves::MutateResidueOP mutres3( new protocols::simple_moves::MutateResidue( 3, "TRP:N_Methylation" ) );
		mutres3->set_update_polymer_dependent( true );
		mutres3->apply(*initial_pose_2chain);
		protocols::simple_moves::MutateResidueOP mutres17( new protocols::simple_moves::MutateResidue( 17, "DTRP:N_Methylation" ) );
		mutres17->set_update_polymer_dependent( true );
		mutres17->apply(*initial_pose_2chain);

		poses_2chain_.clear();
		mirror_poses_2chain_.clear();

		poses_2chain_.push_back(initial_pose_2chain);
		mirror_poses_2chain_.push_back( mirror_pose_with_disulfides( poses_2chain_[1] ) );
		for ( core::Size i=1; i<=23; ++i ) {
			poses_2chain_.push_back( permute_with_disulfides( poses_2chain_[i] ) );
			mirror_poses_2chain_.push_back( mirror_pose_with_disulfides( poses_2chain_[i+1] ) );

			poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(1);
			poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(13);
			poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(12);
			poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(24);
			mirror_poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(1);
			mirror_poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(13);
			mirror_poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(12);
			mirror_poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(24);

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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_, false);
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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the lk_ball_wtd scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_lk_ball_wtd() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::lk_ball_wtd, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing lk_ball_wtd score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_elec scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_elec() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_elec, 1.0 );
		TR << "Testing fa_elec score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_dun scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_dun() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_dun, 1.0 );
		TR << "Testing fa_dun score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
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
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the yhh_planarity scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_yhh_planarity() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::yhh_planarity, 1.0 );
		TR << "Testing yhh_planarity score term." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the full beta_nov15 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_beta_nov15() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->add_weights_from_file("beta_nov15.wts");
		TR << "Testing full beta_nov15 score function." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the full scorefunction that's the current default (whatever that is).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_current_default_scorefxn() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
		TR << "Testing full default score function." << std::endl;
		CyclicGeometryTestHelper helper;
		helper.cyclic_pose_test(scorefxn, poses_2chain_, mirror_poses_2chain_);
		return;
	}

private:
	utility::vector1 < core::pose::PoseOP > poses_2chain_;
	utility::vector1 < core::pose::PoseOP > mirror_poses_2chain_;


};
