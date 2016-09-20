// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file test/core/scoring/SymmDL.cxxtest.hh
/// @brief Unit tests for symmetric pose scoring with the --symmetric_gly_tables option.
/// @detials Mirror image poses (with D- and L-amino acids swapped) should score identically with this option.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pdb1rpb.hh>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

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
#include <protocols/simple_moves/MutateResidue.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/chemical/AA.hh>

//Minimizer
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/MoveMap.hh>

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;


static THREAD_LOCAL basic::Tracer TR("core.scoring.SymmDL_betanov15_Tests.cxxtest");

class SymmDL_betanov15_Tests : public CxxTest::TestSuite {

public:

	void setUp() {
		//core_init();
		core_init_with_additional_options( "-beta_nov15 -score:weights beta_nov15.wts -symmetric_gly_tables true -write_all_connect_info -connect_info_cutoff 0.0" );
	}

	void tearDown() {
	}

	/// @brief Given a residue type, get its mirror-image type.
	///
	core::chemical::ResidueType const & get_mirror_type( core::chemical::ResidueType const &master_type ) {
		if ( !master_type.is_l_aa() && !master_type.is_d_aa() ) return master_type;
		core::chemical::ResidueTypeSet const & residue_type_set = *master_type.residue_type_set();
		if ( master_type.is_l_aa() ) {
			return residue_type_set.name_map( "D"+master_type.name() );
		}
		return residue_type_set.name_map( master_type.name().substr(1) );
	}

	/// @brief Given a pose, flip the L-residues to D-residues.
	///
	void flip_residues( core::pose::Pose &pose ) {
		// Create the new residue and replace it
		for ( core::Size ir=1, irmax=pose.size(); ir<=irmax; ++ir ) {
			core::conformation::ResidueOP new_res = core::conformation::ResidueFactory::create_residue( get_mirror_type( pose.residue(ir).type() ), pose.residue( ir ), pose.conformation());
			core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( ir ), *new_res, pose.conformation(), true );
			pose.replace_residue( ir, *new_res, false );
		}
	}

	/// @brief Given a pose, construct its mirror image.
	///
	void mirror_pose( core::pose::Pose const &master, core::pose::Pose &mirror) {

		for ( core::Size ir=1, irmax=mirror.size(); ir<=irmax; ++ir ) {
			for ( core::Size ia=1, iamax=mirror.residue(ir).type().natoms(); ia<=iamax; ++ia ) {
				core::id::AtomID const curatom(ia, ir);
				numeric::xyzVector<core::Real> xyztemp( master.xyz( curatom ) );
				xyztemp.z( -1.0*xyztemp.z() );
				mirror.set_xyz( curatom, xyztemp );
			}
		}

		mirror.update_residue_neighbors();

		return;
	}

	/// @brief Construct a few L-amino acid poses, and confirm that mirror-image conformations
	/// score identically with a given scorefunction.
	void mirror_pose_test( core::scoring::ScoreFunctionOP sfxn ) {
		core::pose::Pose pose = pdb1rpb_pose();
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

		//Mutate residue 4 to proline.  Since residue 3 is a gly, this makes a nice gly-pro.
		protocols::simple_moves::MutateResidue mutres;
		mutres.set_res_name("PRO");
		mutres.set_target(4);
		mutres.apply(pose);

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
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::cart_bonded, 1.0 );
		TR << "Testing cart_bonded score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_atr scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_fa_atr() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_atr, 1.0 );
		TR << "Testing fa_atr score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_fa_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_rep, 1.0 );
		TR << "Testing fa_rep score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_intra_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_fa_intra_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_rep, 1.0 );
		TR << "Testing fa_intra_rep score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_sol scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_fa_sol() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_sol score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_intra_sol_xover4 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_fa_intra_sol_xover4() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_sol_xover4, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_intra_sol_xover4 score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the lk_ball_wtd scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_lk_ball_wtd() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::lk_ball_wtd, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing lk_ball_wtd score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the fa_elec scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_fa_elec() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_elec, 1.0 );
		TR << "Testing fa_elec score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the pro_close scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_pro_close() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::pro_close, 1.0 );
		TR << "Testing pro_close score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the dslf_fa13 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_dslf_fa13() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::dslf_fa13, 1.0 );
		TR << "Testing dslf_fa13 score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the hbonds scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_hbonds() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
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
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_dun, 1.0 );
		TR << "Testing fa_dun score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the omega scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_omega() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::omega, 1.0 );
		TR << "Testing omega score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the rama scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_rama() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama, 1.0 );
		TR << "Testing rama score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the rama_prepro scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_rama_prepro() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama_prepro, 1.0 );
		TR << "Testing rama_prepro score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the rama_prepro scorefunction and a gly-pro pair.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_rama_prepro2() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama_prepro, 1.0 );
		TR << "Testing rama_prepro score term with a gly-pro pair." << std::endl;
		mirror_pose_test2(scorefxn);
		return;
	}


	/// @brief Tests symmetric scoring with the p_aa_pp scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_p_aa_pp() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::p_aa_pp, 1.0 );
		TR << "Testing p_aa_pp score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the yhh_planarity scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_yhh_planarity() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::yhh_planarity, 1.0 );
		TR << "Testing yhh_planarity score term." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the full beta_nov15 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_beta_nov15() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->add_weights_from_file("beta_nov15.wts");
		TR << "Testing full beta_nov15 score function." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

	/// @brief Tests symmetric scoring with the full default scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_DL_default_scorefxn() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
		TR << "Testing full default score function." << std::endl;
		mirror_pose_test(scorefxn);
		return;
	}

};
