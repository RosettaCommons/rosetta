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

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;


static THREAD_LOCAL basic::Tracer TR("core.scoring.CyclicGeometry_betanov15_Tests.cxxtest");

class CyclicGeometry_betanov15_Tests : public CxxTest::TestSuite {

public:

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
		for ( core::Size ir=1, irmax=pose.n_residue(); ir<=irmax; ++ir ) {
			core::conformation::ResidueOP new_res = core::conformation::ResidueFactory::create_residue( get_mirror_type( pose.residue(ir).type() ), pose.residue( ir ), pose.conformation());
			core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( ir ), *new_res, pose.conformation(), true );
			pose.replace_residue( ir, *new_res, false );
		}
	}

	/// @brief Given a pose, construct its mirror image.
	///
	core::pose::PoseOP
	mirror_pose( core::pose::PoseCOP master) {
		core::pose::PoseOP output_pose( master->clone() );
		flip_residues(*output_pose);

		for ( core::Size ir=1, irmax=output_pose->n_residue(); ir<=irmax; ++ir ) {
			for ( core::Size ia=1, iamax=output_pose->residue(ir).type().natoms(); ia<=iamax; ++ia ) {
				core::id::AtomID const curatom(ia, ir);
				numeric::xyzVector<core::Real> xyztemp( master->xyz( curatom ) );
				xyztemp.z( -1.0*xyztemp.z() );
				output_pose->set_xyz( curatom, xyztemp );
			}
		}

		return output_pose;
	}

	/// @brief Given an input cyclic pose, circularly permute by 1 residue.
	core::pose::PoseOP permute(
		core::pose::PoseOP input_pose
	) {
		core::pose::PoseOP new_pose( new core::pose::Pose );

		for ( core::Size i=2, imax=input_pose->n_residue(); i<=imax; ++i ) {
			core::conformation::ResidueOP new_rsd( input_pose->residue(i).clone() );
			if ( i == 2 ) {
				new_pose->append_residue_by_jump(*new_rsd, 1);
			} else {
				new_pose->append_residue_by_bond(*new_rsd, false);
			}
		}
		new_pose->append_residue_by_bond( *(input_pose->residue(1).clone()), false);
		new_pose->conformation().declare_chemical_bond(1, "N", 9, "C");
		return new_pose;
	}

	void setUp() {
		core_init_with_additional_options( "-beta_nov15 -score:weights beta_nov15.wts -symmetric_gly_tables true -write_all_connect_info -connect_info_cutoff 0.0" );

		// Pull in the cyclic peptide pose (9 residues):
		core::pose::PoseOP initial_pose( new core::pose::Pose );
		core::import_pose::pose_from_file( *initial_pose, "core/scoring/cyclic_peptide.pdb" , core::import_pose::PDB_file );

		// Strip off termini and connect the ends:
		core::pose::remove_variant_type_from_pose_residue( *initial_pose, core::chemical::LOWER_TERMINUS_VARIANT, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose, core::chemical::CUTPOINT_LOWER, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose, core::chemical::CUTPOINT_UPPER, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose, core::chemical::UPPER_TERMINUS_VARIANT, 9 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose, core::chemical::CUTPOINT_LOWER, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose, core::chemical::CUTPOINT_UPPER, 1 );
		initial_pose->conformation().declare_chemical_bond(1, "N", 9, "C");
		initial_pose->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(1);
		initial_pose->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(9);

		poses_.push_back(initial_pose);
		mirror_poses_.push_back( mirror_pose( poses_[1] ) );
		for ( core::Size i=1; i<=8; ++i ) {
			poses_.push_back( permute( poses_[i] ) );
			mirror_poses_.push_back( mirror_pose( poses_[i+1] ) );
		}
	}

	void tearDown() {
	}

	/// @brief Test that the same score is returned for all cyclic permutations.
	///
	void cyclic_pose_test( core::scoring::ScoreFunctionOP sfxn ) {
		//Score all of the poses and confirm that they're all equal to the first
		for ( core::Size i=1; i<=9; ++i ) {
			(*sfxn)(*poses_[i]);
			if ( i>1 ) {
				TR << "\tTesting scoring with circular permutation of " << i - 1 << " residues." << std::endl;
				TS_ASSERT_DELTA(poses_[1]->energies().total_energy(), poses_[i]->energies().total_energy(), std::abs( std::max(poses_[1]->energies().total_energy(), poses_[i]->energies().total_energy())/10000.0 ) );
			}
			//Check mirrored geometry, too:
			TR << "\tTesting scoring with circular permutation of " << i - 1 << " residues and mirroring." << std::endl;
			(*sfxn)(*mirror_poses_[i]);
			TS_ASSERT_DELTA(poses_[1]->energies().total_energy(), mirror_poses_[i]->energies().total_energy(), std::abs( std::max(poses_[1]->energies().total_energy(), mirror_poses_[i]->energies().total_energy())/10000.0 ) );
			TR.flush();
		}

		//Delete the following:
		//poses_[1]->dump_scored_pdb("vcyclic1.pdb", *sfxn);
		//poses_[2]->dump_scored_pdb("vcyclic2.pdb", *sfxn);
		//poses_[3]->dump_scored_pdb("vcyclic3.pdb", *sfxn);
		//poses_[4]->dump_scored_pdb("vcyclic4.pdb", *sfxn);
		//poses_[5]->dump_scored_pdb("vcyclic5.pdb", *sfxn);
		//poses_[6]->dump_scored_pdb("vcyclic6.pdb", *sfxn);
		//poses_[7]->dump_scored_pdb("vcyclic7.pdb", *sfxn);
		//poses_[8]->dump_scored_pdb("vcyclic8.pdb", *sfxn);
		//poses_[9]->dump_scored_pdb("vcyclic9.pdb", *sfxn);

	}

	/// @brief Tests cyclic permutation scoring with the fa_atr scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_atr() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_atr, 1.0 );
		TR << "Testing fa_atr score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_rep, 1.0 );
		TR << "Testing fa_rep score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_intra_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_intra_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_rep, 1.0 );
		TR << "Testing fa_intra_rep score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_sol scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_sol() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_sol score term." << std::endl;
		cyclic_pose_test(scorefxn);
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
		cyclic_pose_test(scorefxn);
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
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_elec scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_elec() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_elec, 1.0 );
		TR << "Testing fa_elec score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the pro_close scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_pro_close() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::pro_close, 1.0 );
		TR << "Testing pro_close score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the dslf_fa13 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_dslf_fa13() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::dslf_fa13, 1.0 );
		TR << "Testing dslf_fa13 score term." << std::endl;
		cyclic_pose_test(scorefxn);
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
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_dun scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_dun() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_dun, 1.0 );
		TR << "Testing fa_dun score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the omega scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_omega() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::omega, 1.0 );
		TR << "Testing omega score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the rama scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_rama() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama, 1.0 );
		TR << "Testing rama score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the rama_prepro scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_rama_prepro() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama_prepro, 1.0 );
		TR << "Testing rama_prepro score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the p_aa_pp scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_p_aa_pp() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::p_aa_pp, 1.0 );
		TR << "Testing p_aa_pp score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the yhh_planarity scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_yhh_planarity() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::yhh_planarity, 1.0 );
		TR << "Testing yhh_planarity score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the full beta_nov15 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_beta_nov15() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->add_weights_from_file("beta_nov15.wts");
		TR << "Testing full beta_nov15 score function." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the full scorefunction that's the current default (whatever that is).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_current_default_scorefxn() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
		TR << "Testing full default score function." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

private:
	utility::vector1 < core::pose::PoseOP > poses_;
	utility::vector1 < core::pose::PoseOP > mirror_poses_;


};
