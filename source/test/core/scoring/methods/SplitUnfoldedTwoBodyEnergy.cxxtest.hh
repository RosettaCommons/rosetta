// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/SplitUnfoldedReferenceEnergy.cxxtest.hh
/// @brief  test suite for SplitUnfoldedReferenceEnergy energy method
/// @author Ron Jacak

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/SplitUnfoldedTwoBodyEnergy.hh>

#include <platform/types.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>


#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoringManager.fwd.hh>
#include <utility/vector1.hh>



// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class SplitUnfoldedReferenceEnergyTests : public CxxTest::TestSuite {

	public:

	PoseOP the_pose;
	Real TOLERATED_ERROR;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		TOLERATED_ERROR = 0.001;
		core_init();

		the_pose = PoseOP( new Pose );
		core::chemical::ResidueTypeSetCOP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		core::pose::make_pose_from_sequence( *the_pose, "DFGLK", *rsd_set );

	}

	// Shared finalization goes here.
	void tearDown() {
		the_pose.reset();
	}


	// --------------- Test Cases --------------- //

	// Test the three different ways of using the SplitUnfoldedReferenceEnergy energy method class. The difference lies in how the class
	// is constructed: with the default constructor, with an empty EMap, or with a full EMap.

	void test_residue_energy_default_file_weights() {
		// AA fa_atr fa_rep fa_sol fa_intra_rep pro_close fa_pair hbond_sr_bb hbond_lr_bb hbond_bb_sc hbond_sc rama omega fa_dun p_aa_pp
		// WEIGHT 0.8 0.44 0.65 0.004 1.0 0.49  0.585 1.17 1.17 1.1 0.2 0.5 0.56 0.32
		// PHE 3.6792 -0.7291 -1.7340 -6.6803 -0.0053 0.0000 0.4257 0.0216 0.0055 0.0000 0.1227 -0.2181 -1.0325 0.1014
		// GLY 1.6906 -0.4879 -1.2009 -0.0181 -0.0033 0.0000 0.3093 0.0139 0.0093 0.0000 -0.6594 -0.1474 0.0000 1.0637

		SplitUnfoldedTwoBodyEnergyOP unweighted_unf_energy( new SplitUnfoldedTwoBodyEnergy( scoring::SPLIT_UNFOLDED_MM, scoring::SPLIT_UNFOLDED_MEDIAN, scoring::UNFOLDED_SPLIT_MM_STD ) );

		EnergyMap emap;
		emap.zero();

		unweighted_unf_energy->residue_energy( the_pose->residue(2), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ split_unfolded_two_body ], -2.4221, TOLERATED_ERROR );

		emap.zero();
		unweighted_unf_energy->residue_energy( the_pose->residue(3), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ split_unfolded_two_body ], 0.0733, TOLERATED_ERROR );

		unweighted_unf_energy = 0;

	}

	void test_residue_energy_no_weights() {

		EnergyMap emap; // same emap used for initalization and holding scores
		emap.zero();

		SplitUnfoldedTwoBodyEnergyOP unweighted_unf_energy( new SplitUnfoldedTwoBodyEnergy( scoring::SPLIT_UNFOLDED_MM, scoring::SPLIT_UNFOLDED_MEDIAN, scoring::UNFOLDED_SPLIT_MM_STD, emap ) );

		unweighted_unf_energy->residue_energy( the_pose->residue(2), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ split_unfolded_two_body ], 0.0, TOLERATED_ERROR );

		emap.zero();
		unweighted_unf_energy->residue_energy( the_pose->residue(3), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ split_unfolded_two_body ], 0.0, TOLERATED_ERROR );

		unweighted_unf_energy = 0;
	}

	void test_residue_energy_emap_weights() {

		EnergyMap emap;

		emap[ fa_atr ] = 0.8;
		emap[ fa_rep ] = 0.634454;
		emap[ fa_sol ] = 1.16497;
		emap[ mm_lj_intra_rep ] = 0.324341;
		emap[ mm_lj_intra_atr ] = 0.537815;
		emap[ mm_twist ] = 0.2662;
		emap[ pro_close ] = 0.374325;
		emap[ fa_elec ] = 0.73;
		emap[ hbond_sr_bb ] = 0.656728;
		emap[ hbond_lr_bb ] = 1.50186;
		emap[ hbond_bb_sc ] = 1.45367;
		emap[ hbond_sc ] = 1.18477;
		emap[ dslf_fa13 ] = 1;
		emap[ unfolded ] = -0.4;
		emap[ split_unfolded_two_body ] = 0.6;

		SplitUnfoldedTwoBodyEnergyOP weighted_unf_energy_emap( new SplitUnfoldedTwoBodyEnergy( scoring::SPLIT_UNFOLDED_MM, scoring::SPLIT_UNFOLDED_MEDIAN, scoring::UNFOLDED_SPLIT_MM_STD, emap ) );

		emap.zero();
		weighted_unf_energy_emap->residue_energy( the_pose->residue(4), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ split_unfolded_two_body ], -1.4525, TOLERATED_ERROR );

		emap.zero();
		weighted_unf_energy_emap->residue_energy( the_pose->residue(5), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ split_unfolded_two_body ], 0.2137, TOLERATED_ERROR );

		emap.zero();
		weighted_unf_energy_emap->residue_energy( the_pose->residue(1), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ split_unfolded_two_body ], 0.7131, TOLERATED_ERROR );

		weighted_unf_energy_emap = 0;
	}

};
