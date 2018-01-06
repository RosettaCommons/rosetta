// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/ReadPoseExtraScoreFilter.cxxtest.hh
/// @brief  test for ReadPoseExtraScoreFilter
/// @author Jack Maguire, jack@med.unc.edu

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
//#include <test/util/rosettascripts.hh>
//#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/simple_filters/ReadPoseExtraScoreFilter.hh>
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>

// Utility Headers
//#include <utility/file/FileName.hh>
//#include <basic/Tracer.hh>

//static basic::Tracer TR("protocols.simple_filters.ReadPoseExtraScoreFilter.cxxtest.hh");

// --------------- Test Class --------------- //

class ReadPoseExtraScoreFilter : public CxxTest::TestSuite {

public:

	void setUp() {
	}

	void tearDown() {
	}

	void test_use_filter() {
		using namespace core::pose;

		Pose crash_test_dummy;

		//extra_pose_info_util takes string addresses, so we have to declare them
		std::string const name1( "testA" );
		std::string const name2( "testB" );
		std::string const name3( "testC" );
		std::string const value3( "1.5" );

		setPoseExtraScore( crash_test_dummy, name1, -1.0 );
		setPoseExtraScore( crash_test_dummy, name2, 0.0 );
		setPoseExtraScore( crash_test_dummy, name3, value3 );

		protocols::simple_filters::ReadPoseExtraScoreFilter filter;
		filter.set_threshold( 0.5 );
		filter.set_term_name( name1 );
		TS_ASSERT_EQUALS( filter.compute( crash_test_dummy ), -1.0 );
		TS_ASSERT( filter.apply( crash_test_dummy ) );

		filter.set_term_name( name2 );
		TS_ASSERT_EQUALS( filter.compute( crash_test_dummy ), 0.0 );
		TS_ASSERT( filter.apply( crash_test_dummy ) );

		filter.set_term_name( name3 );
		TS_ASSERT_EQUALS( filter.compute( crash_test_dummy ), 1.5 );
		TS_ASSERT( ! filter.apply( crash_test_dummy ) );

	}

};
