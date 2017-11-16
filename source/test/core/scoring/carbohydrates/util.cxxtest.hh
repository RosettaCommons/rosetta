// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/pose/carbohydrates/util.cxxtest.hh
/// @brief   Test suite for utility functions for scoring carbohydrate-containing poses
/// @author  Labonte <JWLabonte@jhu.edu>


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Header
#include <core/scoring/carbohydrates/util.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/types.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Basic Headers
#include <basic/Tracer.hh>


static basic::Tracer TR( "core.scoring.carbohydrates.util.cxxtest" );


class CarbohydrateScoringUtilityFunctionTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		using namespace core::pose;

		core_init_with_additional_options( "-include_sugars" );

		// This is a completely unnatural sugar, created solely to stretch and test the system.
		make_pose_from_saccharide_sequence( test_sugar_, "b-D-Glcp-(1->8)-[a-D-Glcp-(1->4)]-a-Neup-(2->"
			"4)-[a-D-Glcp-(1->6)]-b-D-Fruf-(2->6)-a-D-Galp-(1->6)-b-D-Manp-(1->"
			"4)-[a-D-Glcp-(1->3)]-b-L-Gulp-(1->4)-b-D-Galp-(1->3)-a-D-Glcp-(1->4)-D-Glcp" );
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	// Confirm that the helper functions to select the parameters for the various functional forms for scoring
	// glycosidic torsions work correctly.
	void test_function_assignment_based_on_type_of_linakge()
	{
		using namespace core::id;
		using namespace core::scoring::carbohydrates;

		// The test sugar main chain:
		// b-D-Glcp-(1->8)-a-Neup-(2->4)-b-D-Fruf-(2->6)-a-D-Galp-(1->6)-b-D-Manp-(1->4)-b-L-Gulp-(1->
		// 4)-b-D-Galp-(1->3)-a-D-Glcp-(1->4)-D-Glcp

		// Torsions tested from left to right:

		// b-D-Glcp-(1->8)-Neup (Fruf is not a pyranose.)
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			phi_dihedral, test_sugar_, 9 ), BETA_LINKS );
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			psi_dihedral, test_sugar_, 9 ), BETA6_LINKS );  // (really, for beta-anything-exocyclic)
		TS_ASSERT_EQUALS( get_omega_preference_for_residue_in_pose( test_sugar_, 9 ), PREFERENCE_NA );

		// a-Neup-(2->4)-D-Fruf (Fruf is not a pyranose.)
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			phi_dihedral, test_sugar_, 8 ), ALPHA_LINKS );
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			psi_dihedral, test_sugar_, 8 ), LINKAGE_NA );
		TS_ASSERT_EQUALS( get_omega_preference_for_residue_in_pose( test_sugar_, 8 ), PREFERENCE_NA );

		// b-D-Fruf-(2->6)-D-Galp (Fruf is not a pyranose; D-Galp's O4 is axial in 4C1.)
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			phi_dihedral, test_sugar_, 7 ), LINKAGE_NA );
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			psi_dihedral, test_sugar_, 7 ), LINKAGE_NA );
		TS_ASSERT_EQUALS( get_omega_preference_for_residue_in_pose( test_sugar_, 7 ), ANTI );

		// a-D-Galp-(1->6)-D-Manp (D-Manp's O4 is equatorial in 4C1.)
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			phi_dihedral, test_sugar_, 6 ), ALPHA_LINKS );
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			psi_dihedral, test_sugar_, 6 ), ALPHA6_LINKS );
		TS_ASSERT_EQUALS( get_omega_preference_for_residue_in_pose( test_sugar_, 6 ), GAUCHE_EFFECT );

		// b-D-Manp-(1->4)-L-Gulp (L-Gulp's O4 is axial in 1C4.)
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			phi_dihedral, test_sugar_, 5 ), BETA_LINKS );
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			psi_dihedral, test_sugar_, 5 ), _2AX_3EQ_4AX_LINKS );
		TS_ASSERT_EQUALS( get_omega_preference_for_residue_in_pose( test_sugar_, 5 ), PREFERENCE_NA );

		// b-L-Gulp-(1->4)-D-Galp (D-Galp's O4 is axial in 4C1.)
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			phi_dihedral, test_sugar_, 4 ), BETA_LINKS );
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			psi_dihedral, test_sugar_, 4 ), _2AX_3EQ_4AX_LINKS );
		TS_ASSERT_EQUALS( get_omega_preference_for_residue_in_pose( test_sugar_, 4 ), PREFERENCE_NA );

		// b-D-Galp-(1->3)-D-Glcp (All of D-Glcp's substituents are equatorial in 4C1.)
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			phi_dihedral, test_sugar_, 3 ), BETA_LINKS );
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			psi_dihedral, test_sugar_, 3 ), _2AX_3EQ_4AX_LINKS );
		TS_ASSERT_EQUALS( get_omega_preference_for_residue_in_pose( test_sugar_, 3 ), PREFERENCE_NA );

		// a-D-Glcp-(1->4)-D-Glcp (All of D-Glcp's substituents are equatorial in 4C1.)
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			phi_dihedral, test_sugar_, 2 ), ALPHA_LINKS );
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			psi_dihedral, test_sugar_, 2 ), _2EQ_3AX_4EQ_LINKS );
		TS_ASSERT_EQUALS( get_omega_preference_for_residue_in_pose( test_sugar_, 2 ), PREFERENCE_NA );

		// The reducing end has no phi, psi, or omega.
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			phi_dihedral, test_sugar_, 1 ), LINKAGE_NA );
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			psi_dihedral, test_sugar_, 1 ), LINKAGE_NA );
		TS_ASSERT_EQUALS( get_omega_preference_for_residue_in_pose( test_sugar_, 1 ), PREFERENCE_NA );

		// The test sugar branches:
		// Branch from Gul4: [a-D-Glcp-(1->3)]-

		// a-D-Glcp-(1->3)-L-Gulp (L-Gulp's O3 is axial in 1C4.)
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			phi_dihedral, test_sugar_, 10 ), ALPHA_LINKS );
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			psi_dihedral, test_sugar_, 10 ), _2EQ_3AX_4EQ_LINKS );
		TS_ASSERT_EQUALS( get_omega_preference_for_residue_in_pose( test_sugar_, 10 ), PREFERENCE_NA );

		// Branch from Fru7: [a-D-Glcp-(1->6)]-

		// a-D-Glcp-(1->6)-D-Fruf (Fruf is not a pyranose.)
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			phi_dihedral, test_sugar_, 11 ), ALPHA_LINKS );
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			psi_dihedral, test_sugar_, 11 ), ALPHA6_LINKS );
		TS_ASSERT_EQUALS( get_omega_preference_for_residue_in_pose( test_sugar_, 11 ), PREFERENCE_NA );

		// Branch from Neup8: [a-D-Glcp-(1->4)]-

		// a-D-Glcp-(1->4)-Neup (Neup is a 2-ketose, so its numbers need to be shifted for psi linkage types.)
		TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
			phi_dihedral, test_sugar_, 12 ), ALPHA_LINKS );
		//TS_ASSERT_EQUALS( get_CHI_energy_function_linkage_type_for_residue_in_pose(
		//  psi_dihedral, test_sugar_, 12 ), _2AX_3EQ_4AX_LINKS );  // TODO: FIX
		TS_ASSERT_EQUALS( get_omega_preference_for_residue_in_pose( test_sugar_, 12 ), PREFERENCE_NA );
	}


private:  // Private datum ////////////////////////////////////////////////////
	core::pose::Pose test_sugar_;
};  // class CarbohydrateScoringUtilityFunctionTests
