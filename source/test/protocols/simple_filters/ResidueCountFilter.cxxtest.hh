// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/ResidueCountFilter.cxxtest.hh
/// @brief  test for ResidueCountFilter mover
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <protocols/simple_filters/ResidueCountFilter.hh>

// Platform Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>


// Utility Headers
#include <basic/Tracer.hh>

// C++ Headers
#include <sstream>

static basic::Tracer TR("protocols.simple_filters.ResidueCountFilter.cxxtest.hh");

// --------------- Test Class --------------- //

class ResidueCountFilterTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();


	}

	void tearDown() {
	}

	void test_residue_count_filter_manual() {
		using namespace protocols::simple_filters;

		core::pose::Pose pose(create_twores_1ubq_pose());

		ResidueCountFilter rcf;

		TS_ASSERT(rcf.compute(pose) == 2);

		// max
		rcf.enable_max_residue_count(true);
		rcf.max_residue_count(1);
		TS_ASSERT(rcf.max_residue_count() == 1);
		TS_ASSERT(rcf.enable_max_residue_count() == true);
		TS_ASSERT(!rcf.apply(pose));

		rcf.max_residue_count(2);
		TS_ASSERT(rcf.apply(pose));

		rcf.max_residue_count(3);
		TS_ASSERT(rcf.apply(pose));

		rcf.enable_max_residue_count(false);
		rcf.max_residue_count(0);
		TS_ASSERT(rcf.apply(pose));


		// min
		rcf.enable_min_residue_count(true);
		rcf.min_residue_count(3);
		TS_ASSERT(rcf.min_residue_count() == 3);
		TS_ASSERT(rcf.enable_min_residue_count() == true);
		TS_ASSERT(!rcf.apply(pose));

		rcf.min_residue_count(2);
		TS_ASSERT(rcf.apply(pose));

		rcf.min_residue_count(1);
		TS_ASSERT(rcf.apply(pose));

		rcf.enable_min_residue_count(false);
		rcf.min_residue_count(4);
		TS_ASSERT(rcf.apply(pose));
	}

	void test_residue_count_report() {

		using namespace protocols::simple_filters;

		core::pose::Pose pose(create_twores_1ubq_pose());

		ResidueCountFilter rcf;

		std::stringstream str_stream;
		rcf.report(str_stream, pose);
		TS_ASSERT_EQUALS( str_stream.str(), "Residue Count: 2\n");

		TS_ASSERT_EQUALS(rcf.report_sm(pose), 2.0);
	}
};
