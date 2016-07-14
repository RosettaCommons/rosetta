// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/constraint_generator/ConstraintsManagerTests.cxxtest.hh
/// @brief  unit tests for constraint generator caching functions
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <protocols/constraint_generator/ConstraintsManager.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

// Boost headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR("ConstraintsManagerTests");


class ConstraintsManagerTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();
	}

	void tearDown(){

	}

	void test_pose_cache()
	{
		using core::scoring::constraints::ConstraintOP;
		using core::scoring::constraints::ConstraintCOPs;

		protocols::constraint_generator::ConstraintsManager const & manager =
			*protocols::constraint_generator::ConstraintsManager::get_instance();

		core::pose::Pose pose = create_trpcage_ideal_pose();

		std::string const csts_name = "TestNullConstraints";

		// first, try to get things from pose. this should throw an error because there is no map
		TS_ASSERT( !manager.has_stored_constraints( pose, csts_name ) );
		TS_ASSERT_THROWS_ANYTHING( manager.retrieve_constraints( pose, csts_name ) );

		// This should setup blank data in the pose and not throw an error
		ConstraintCOPs csts;

		// create silly constraint list with two null constraints
		csts.clear();
		csts.push_back( ConstraintOP() );
		csts.push_back( ConstraintOP() );

		manager.store_constraints( pose, csts_name, csts );

		TS_ASSERT( manager.has_stored_constraints( pose, csts_name ) );
		ConstraintCOPs test_csts = manager.retrieve_constraints( pose, csts_name );
		TS_ASSERT_EQUALS( test_csts, csts );

		// should be able to add second set, and not disrupt the first
		std::string const csts_name2 = "TestNullConstraints2";

		ConstraintCOPs const csts2 = boost::assign::list_of(ConstraintOP())(ConstraintOP())(ConstraintOP());

		// first, try to get things from pose. This should throw an error because there are no csts with that name
		TS_ASSERT( !manager.has_stored_constraints( pose, csts_name2 ) );
		TS_ASSERT_THROWS_ANYTHING( manager.retrieve_constraints( pose, csts_name2 ) );

		manager.store_constraints( pose, csts_name2, csts2 );
		TS_ASSERT_EQUALS( manager.retrieve_constraints( pose, csts_name2 ), csts2 );
		TS_ASSERT_EQUALS( manager.retrieve_constraints( pose, csts_name ), csts );
		TS_ASSERT_DIFFERS( manager.retrieve_constraints( pose, csts_name ), csts2 );

		// now, remove first set, second set should still be there
		manager.remove_constraints( pose, csts_name );
		TS_ASSERT( !manager.has_stored_constraints( pose, csts_name ) );
		TS_ASSERT_EQUALS( manager.retrieve_constraints( pose, csts_name2 ), csts2 );

		// remove second set, and everything should be empty
		manager.remove_constraints( pose, csts_name2 );
		TS_ASSERT( !manager.has_stored_constraints( pose, csts_name ) );
		TS_ASSERT( !manager.has_stored_constraints( pose, csts_name2 ) );
	}

};

