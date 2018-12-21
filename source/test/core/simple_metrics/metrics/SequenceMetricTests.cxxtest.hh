// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/simple_metrics/metrics/SequenceMetricTests.cxxtest.hh
/// @brief  Unit tests for the SequenceMetric simple metric.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/simple_metrics/metrics/SequenceMetric.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("SequenceMetricTests");


class SequenceMetricTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

		core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );
		core::pose::make_pose_from_sequence( *pose, "RSTLNEX[ORN]X[DORN]Y[DTYR]YS[DSER:CtermProteinFull]", core::chemical::FA_STANDARD, true, false );
		pose_ = pose;
	}

	void tearDown(){}



	void test_oneletter_code(){
		using namespace core::simple_metrics::metrics;
		SequenceMetric metric1, metric2;
		metric1.set_output_mode( "oneletter" );
		metric2.set_output_mode( SMM_ONELETTER_CODE );

		std::string const expected_result( "RSTLNEXXYYS" );

		TS_ASSERT_EQUALS( metric1.calculate( *pose_ ), expected_result);
		TS_ASSERT_EQUALS( metric2.calculate( *pose_ ), expected_result);
	}

	void test_threeletter_code(){
		using namespace core::simple_metrics::metrics;
		SequenceMetric metric1, metric2;
		metric1.set_output_mode( "threeletter" );
		metric2.set_output_mode( SMM_THREELETTER_CODE );

		std::string const expected_result( "ARG,SER,THR,LEU,ASN,GLU,ORN,ORN,DTY,TYR,DSE" );

		TS_ASSERT_EQUALS( metric1.calculate( *pose_ ), expected_result);
		TS_ASSERT_EQUALS( metric2.calculate( *pose_ ), expected_result);
	}

	void test_base_name(){
		using namespace core::simple_metrics::metrics;
		SequenceMetric metric1, metric2;
		metric1.set_output_mode( "basename" );
		metric2.set_output_mode( SMM_BASE_NAME );

		std::string const expected_result( "ARG,SER,THR,LEU,ASN,GLU,ORN,DORN,DTYR,TYR,DSER" );

		TS_ASSERT_EQUALS( metric1.calculate( *pose_ ), expected_result);
		TS_ASSERT_EQUALS( metric2.calculate( *pose_ ), expected_result);
	}

	void test_full_name(){
		using namespace core::simple_metrics::metrics;
		SequenceMetric metric1, metric2;
		metric1.set_output_mode( "fullname" );
		metric2.set_output_mode( SMM_FULL_NAME );

		std::string const expected_result( "ARG:NtermProteinFull,SER,THR,LEU,ASN,GLU,ORN,DORN,DTYR,TYR,DSER:CtermProteinFull" );

		TS_ASSERT_EQUALS( metric1.calculate( *pose_ ), expected_result);
		TS_ASSERT_EQUALS( metric2.calculate( *pose_ ), expected_result);
	}

private:

	core::pose::PoseCOP pose_;

};
