// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/io/merge_and_split_behaviors_io.cxxtest.hh
/// @brief   Test suite for residue-splitting-behaviors database loading
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/io/merge_and_split_behaviors_io.hh>

// Basic header
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// C++ header
#include <utility>

static basic::Tracer TR( "core.io.merge_and_split_behaviors_io" );

class MergeAndSplitBehaviorsIOTests : public CxxTest::TestSuite {
public: // Standard methods ///////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init();
	}

	// Destruction
	void tearDown()
	{}


public: // Tests //////////////////////////////////////////////////////////////
	// Confirm that alternative PDB 3-letter-codes are loaded correctly from the database.
	void test_read_merge_behaviors_from_database_file()
	{
		using namespace std;
		using namespace core::chemical::io;

		TR <<  "Testing read_merge_behaviors_from_database_file() method."  << endl;

		MergeBehaviorMap voltron( read_merge_behaviors_from_database_file(
			"core/io/merge_residue_behaviors.fake" ) );

		TS_ASSERT_EQUALS( voltron.size(), 5 );
		ResidueMergeInstructions const & core( voltron[ "Black Lion" ] );
		ResidueMergeInstructions const & right_arm( voltron[ "Red Lion" ] );
		ResidueMergeInstructions const & left_arm( voltron[ "Green Lion" ] );
		ResidueMergeInstructions const & right_leg( voltron[ "Blue Lion" ] );
		ResidueMergeInstructions const & left_leg( voltron[ "Yellow Lion" ] );

		TS_ASSERT_EQUALS( core.first, mrb_merge_w_next );
		TS_ASSERT_EQUALS( core.second.size(), 4 );
		TS_ASSERT_EQUALS( core.second.at( "mouth" ), "face" );
		TS_ASSERT_EQUALS( right_arm.first, mrb_merge_w_prev );
		TS_ASSERT_EQUALS( right_arm.second.at( "mouth" ), "hand" );
		TS_ASSERT_EQUALS( left_arm.first, mrb_merge_w_prev );
		TS_ASSERT_EQUALS( left_arm.second.at( "mouth" ), "hand" );
		TS_ASSERT_EQUALS( right_leg.first, mrb_merge_w_prev );
		TS_ASSERT_EQUALS( right_leg.second.at( "mouth" ), "foot" );
		TS_ASSERT_EQUALS( left_leg.first, mrb_merge_w_prev );
		TS_ASSERT_EQUALS( left_leg.second.at( "mouth" ), "foot" );
	}

	// Confirm that alternative PDB 3-letter-codes are loaded correctly from the database.
	void test_read_split_behaviors_from_database_file()
	{
		using namespace std;
		using namespace core::chemical::io;

		TR <<  "Testing read_split_behaviors_from_database_file() method."  << endl;

		SplitBehaviorsMap behaviors( read_split_behaviors_from_database_file(
			"core/io/split_residue_behaviors.fake" ) );

		TS_ASSERT_EQUALS( behaviors.size(), 1 );
		SplitBehaviors const & behavior( behaviors[ "2FC" ] );

		TS_ASSERT_EQUALS( behavior.first[ 1 ].first, "HD " );
		TS_ASSERT_EQUALS( behavior.first[ 1 ].second, "Harvey Dent" );
		TS_ASSERT_EQUALS( behavior.first[ 2 ].first, "2Fc" );
		TS_ASSERT_EQUALS( behavior.first[ 2 ].second, "Two-Face" );
		TS_ASSERT_EQUALS( behavior.second[ 1 ].size(), 7 );
		TS_ASSERT_EQUALS( behavior.second[ 2 ].size(), 7 );
		TS_ASSERT_EQUALS( behavior.second[ 1 ].at( "RARM" ), " ARM" );
		TS_ASSERT_EQUALS( behavior.second[ 2 ].at( "LARM" ), " ARM" );
		TS_ASSERT_EQUALS( behavior.second[ 1 ].at( "RFT " ), "FOOT" );
	}
};  // class SplitBehaviorsIOTests
