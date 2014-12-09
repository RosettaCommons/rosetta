// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/match/ScaffoldBuildPoint.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>

#include <core/conformation/Residue.hh>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers
#include <utility/exit.hh>

/// Project headers
#include <core/types.hh>

// C++ headers
#include <string>
#include <iostream>

//Auto Headers
#include <utility/vector1.hh>


using namespace protocols::match;
using namespace protocols::match::upstream;


// --------------- Test Class --------------- //

class ScaffoldBuildPointTests : public CxxTest::TestSuite {

	public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.


	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	void test_scaffold_build_point_ctor() {
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		OriginalBackboneBuildPoint res2bp( trpcage.residue( 2 ), 1 );
		TS_ASSERT( res2bp.N_pos().distance_squared(  trpcage.residue( 2 ).xyz( "N"  ) ) < 1e-6 );
		TS_ASSERT( res2bp.CA_pos().distance_squared( trpcage.residue( 2 ).xyz( "CA" ) ) < 1e-6 );
		TS_ASSERT( res2bp.C_pos().distance_squared(  trpcage.residue( 2 ).xyz( "C"  ) ) < 1e-6 );
		TS_ASSERT( res2bp.O_pos().distance_squared(  trpcage.residue( 2 ).xyz( "O"  ) ) < 1e-6 );
	}


};
