// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/simple_metrics/metrics/CustomRealValueMetricTests.cxxtest.hh
/// @brief  Unit tests for the CustomRealValueMetric, which lets you store arbitrary floats in a Pose.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/simple_metrics/metrics/CustomRealValueMetric.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/simple_metrics/util.hh>
#include <core/simple_metrics/SimpleMetricData.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility Headers

// Basic Headers
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>

static basic::Tracer TR("CustomRealValueMetricTests");


class CustomRealValueMetricTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}

	void tearDown() {

	}


	void test_custom_real_value_metric() {
		TR << "Starting unit test CustomRealValueMetricTests::test_custom_real_value_metric." << std::endl;
		core::pose::PoseOP pose( core::import_pose::pose_from_file( "core/simple_metrics/metrics/testpose.pdb" ) );

		core::simple_metrics::metrics::CustomRealValueMetric metric1, metric2, metric3;
		metric1.set_value( 5.3 );
		metric2.set_value( 8.6 );
		metric3.set_value( -9.1 );

		metric1.apply( "first_metric", *pose );
		metric2.apply( "second_metric", *pose );

		TS_ASSERT( pose->data().has( core::pose::datacache::CacheableDataType::SIMPLE_METRIC_DATA ) );

		core::simple_metrics::SimpleMetricDataCOP data( core::simple_metrics::get_sm_data( *pose ) );
		TS_ASSERT( data != nullptr );
		TS_ASSERT( data->get_real_metric_data().count( "first_metric" ) == 1 );
		TS_ASSERT( data->get_real_metric_data().count( "second_metric" ) == 1 );
		TS_ASSERT_DELTA( data->get_real_metric_data().at("first_metric"), 5.3, 0.00001 );
		TS_ASSERT_DELTA( data->get_real_metric_data().at("second_metric"), 8.6, 0.00001 );

		TR << "The following should result in a throw:" << std::endl;
		TS_ASSERT_THROWS_ANYTHING( metric3.apply("first_metric", *pose) );
		TR << "An error message should appear above this line if the unit test is working correctly." << std::endl;

		TR << "Completed unit test CustomRealValueMetricTests::test_custom_real_value_metric." << std::endl;
	}


};
