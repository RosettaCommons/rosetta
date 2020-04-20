// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/filters/FilterValueMetricTest.cxxtest.hh
/// @brief  Test the FilterValueMetric class
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>

// Project Headers
#include <protocols/filters/FilterValueMetric.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("FilterValueMetricTests");


class FilterValueMetricTests : public CxxTest::TestSuite {

public:

	void setUp(){
		core_init();
	}

	void tearDown(){
	}


	void test_constructor() {
		StubMultiFilterOP value_filter = utility::pointer::make_shared< StubMultiFilter >();
		value_filter->set( {43, 12, 42} );

		protocols::filters::FilterValueMetric filter_metric( value_filter );

		core::pose::Pose pose;
		TS_ASSERT_EQUALS( filter_metric.calculate(pose), 43 );
		TS_ASSERT_EQUALS( filter_metric.calculate(pose), 12 );
		TS_ASSERT_EQUALS( filter_metric.calculate(pose), 42 );
	}

};
