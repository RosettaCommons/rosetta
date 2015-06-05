// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/SplitUnfoldedTwoBodyPotential.cxxtest.hh
/// @brief  test suite for core::scoring::SplitUnfoldedTwoBodyPotential.cc
/// @author Ron Jacak

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/SplitUnfoldedTwoBodyPotential.hh>

#include <platform/types.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <core/pose/Pose.hh>


#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>

// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>

// AUTO-REMOVED #include <basic/options/option.hh>

#include <basic/database/open.hh>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>



// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::chemical;

class SplitUnfoldedTwoBodyPotentialTests : public CxxTest::TestSuite {

	public:

	PoseOP pose;
	SplitUnfoldedTwoBodyPotentialOP tcre_potential;
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

		pose = PoseOP( new Pose );
		core::chemical::ResidueTypeSetCOP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		core::pose::make_pose_from_sequence( *pose, "DFGLK", *rsd_set );

		tcre_potential = SplitUnfoldedTwoBodyPotentialOP( new SplitUnfoldedTwoBodyPotential( basic::database::full_name( "scoring/score_functions/split_unfolded/mm_database_median" ) ) );

	}

	// Shared finalization goes here.
	void tearDown() {
		pose.reset();
		tcre_potential.reset();
	}


	// --------------- Test Cases --------------- //
	void test_get_restype_emap() {

		// atr  rep  sol  intra_rep  pro_close  pair  hb_sr_bb  hb_lr_bb  hb_bb_sc  hb_sc  rama  omega  dun  paapp
		// GLY 1.6906 -0.4879 -1.2009 -0.0181 -0.0033 0.0000 0.3093 0.0139 0.0093 0.0000 -0.6594 -0.1474 0.0000 1.0637
		// LEU 3.2405 -0.6543 -1.7746 -2.0520 -0.0049 0.0000 0.5425 0.0343 0.0046 0.0000 0.1254 -0.1717 -0.6123 0.1728
		
		EnergyMap emap;
		emap.zero();
		tcre_potential->get_restype_emap( (pose->residue(3)).type(), emap );
		TS_ASSERT_DELTA( emap[ fa_atr ], -2.8032, TOLERATED_ERROR );
		TS_ASSERT_DELTA( emap[ fa_sol ], 1.8494, TOLERATED_ERROR );
		TS_ASSERT_DELTA( emap[ mm_lj_intra_rep ], 0.0, TOLERATED_ERROR );
		TS_ASSERT_DELTA( emap[ mm_lj_intra_atr ], 0.0, TOLERATED_ERROR );

		emap.zero();
		tcre_potential->get_restype_emap( (pose->residue(4)).type(), emap );
		TS_ASSERT_DELTA( emap[ fa_rep ], 0.3885, TOLERATED_ERROR );
		TS_ASSERT_DELTA( emap[ fa_elec ], -0.2660, TOLERATED_ERROR );

	}

};


