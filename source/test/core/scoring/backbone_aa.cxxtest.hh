// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file test/core/scoring/backbone_aa.cxxtest.hh
/// @brief Unit tests for the BACKBONE_AA lines in params file.
/// @detials The rama, p_aa_pp, and rama_prepro scores for residues with BACKBONE_AA lines should match the corresponding scores for the corresponding residues in the corresponding geometries.
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


static basic::Tracer TR("core.scoring.BackboneAATests.cxxtest");

class BackboneAATests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options( "-symmetric_gly_tables true -write_all_connect_info -connect_info_cutoff 0.0 -extra_res_fa core/scoring/fake_restype_with_tyr_backbone_aa.params" );

		core::pose::PoseOP tyrpose( new core::pose::Pose );
		core::pose::make_pose_from_sequence(*tyrpose, "GYG", "fa_standard", true);
		tyrpose->set_omega(1, 180);
		tyrpose->set_omega(2, 180);

		core::pose::PoseOP ncaapose( tyrpose->clone() );
		core::chemical::ResidueTypeSetCOP standard_residues( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
		core::conformation::ResidueOP new_rsd1( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("FAKE_RESTYPE") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( ncaapose->residue( 2 ), *new_rsd1, ncaapose->conformation(), true);
		ncaapose->replace_residue( 2, *new_rsd1, false );

		tyr_poses_.clear();

		tyrpose->set_phi(2,-61.2);
		tyrpose->set_psi(2,-43.2);
		tyr_poses_.push_back(tyrpose->clone());
		tyrpose->set_phi(2,-133.7);
		tyrpose->set_psi(2,132.3);
		tyr_poses_.push_back(tyrpose->clone());
		tyrpose->set_phi(2,-59.6);
		tyrpose->set_psi(2,-59.9);
		tyr_poses_.push_back(tyrpose->clone());
		tyrpose->set_phi(2,-62.3);
		tyrpose->set_psi(2,-20.1);
		tyr_poses_.push_back(tyrpose->clone());
		tyrpose->set_phi(2,-150.3);
		tyrpose->set_psi(2,-179.3);
		tyr_poses_.push_back(tyrpose->clone());
		tyrpose->set_phi(2,-120.2);
		tyrpose->set_psi(2,-160.3);
		tyr_poses_.push_back(tyrpose->clone());
		tyrpose->set_phi(2,-53.2);
		tyrpose->set_psi(2,139.3);
		tyr_poses_.push_back(tyrpose->clone());
		tyrpose->set_phi(2,-32.1);
		tyrpose->set_psi(2,-160.3);
		tyr_poses_.push_back(tyrpose->clone());
		tyrpose->set_phi(2,-80.7);
		tyrpose->set_psi(2,-70.3);
		tyr_poses_.push_back(tyrpose->clone());
		tyrpose->set_phi(2,-2.1);
		tyrpose->set_psi(2,-89.3);
		tyr_poses_.push_back(tyrpose->clone());

		for ( core::Size i=1; i<=10; ++i ) {
			tyrpose->set_phi( 2, -1.0*tyr_poses_[i]->phi(2) );
			tyrpose->set_psi( 2, -1.0*tyr_poses_[i]->psi(2) );
			tyr_poses_.push_back(tyrpose->clone());
		}

		ncaa_poses_.clear();
		for ( core::Size i=1; i<=20; ++i ) {
			ncaapose->set_phi( 2, tyr_poses_[i]->phi(2) );
			ncaapose->set_psi( 2, tyr_poses_[i]->psi(2) );
			ncaa_poses_.push_back( ncaapose->clone() );
		}

	}

	void tearDown() {
	}

	/// @brief Test that the same score is returned for the BACKBONE_AA restype as would be returned for tyrosine, for a given score term.
	///
	void backbone_aa_test( core::scoring::ScoreFunctionOP sfxn ) {
		//Score all of the poses and confirm that they're all equal to the corresponding poses with tyrosine
		TR << "Scoring " << tyr_poses_.size() << " tyrosine poses and " << ncaa_poses_.size() << " ncaa poses." << std::endl;
		for ( core::Size i=1; i<=20; ++i ) {
			(*sfxn)(*tyr_poses_[i]);
			(*sfxn)(*ncaa_poses_[i]);
			TS_ASSERT_DELTA(tyr_poses_[i]->energies().total_energy(), ncaa_poses_[i]->energies().total_energy(), std::abs( std::max(tyr_poses_[i]->energies().total_energy(), ncaa_poses_[i]->energies().total_energy())/10000.0 ) );
		}

	}

	/// @brief Tests cyclic permutation scoring with the rama scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_backbone_aa_rama() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama, 1.0 );
		TR << "Testing rama score term." << std::endl;
		backbone_aa_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the rama_prepro scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_backbone_aa_rama_prepro() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama_prepro, 1.0 );
		TR << "Testing rama_prepro score term." << std::endl;
		backbone_aa_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the p_aa_pp scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_backbone_aa_p_aa_pp() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::p_aa_pp, 1.0 );
		TR << "Testing p_aa_pp score term." << std::endl;
		backbone_aa_test(scorefxn);
		return;
	}

private:
	utility::vector1 < core::pose::PoseOP > tyr_poses_;
	utility::vector1 < core::pose::PoseOP > ncaa_poses_;

};
