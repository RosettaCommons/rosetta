// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/filters/ReplicateFilter.cxxtest.hh
/// @brief  test for the replicate filter
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/rosettascripts.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/filters/ReplicateFilter.hh>

// Utility Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.filters.ReplicateFilter.cxxtest.hh");

// --------------- Test Class --------------- //

class ReplicateFilterTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP testpose_;
public:

	void setUp() {
		core_init();

		testpose_ = create_twores_1ubq_poseop(); // Identity doesn't matter
	}

	void tearDown() {
	}

	void test_thresholding() {
		protocols::filters::FilterOP sfP( new StubFilter( true, 222) );
		protocols::filters::FilterOP sfN( new StubFilter( true, -335) );

		protocols::filters::ReplicateFilter filterP(sfP,1);
		protocols::filters::ReplicateFilter filterN(sfN,1);
		//default 0 threshold
		TS_ASSERT( ! filterP.apply(*testpose_) );
		TS_ASSERT( filterN.apply(*testpose_) );
		filterP.threshold( 1000 );
		filterN.threshold( 1000 );
		TS_ASSERT( filterP.apply(*testpose_) );
		TS_ASSERT( filterN.apply(*testpose_) );
		filterP.threshold( -1000 );
		filterN.threshold( -1000 );
		TS_ASSERT( ! filterP.apply(*testpose_) );
		TS_ASSERT( ! filterN.apply(*testpose_) );
	}

	void test_basic() {
		utility::vector1<core::Real> values;
		values.push_back(-2);
		values.push_back(-4);
		values.push_back(-1);
		values.push_back(-6);
		values.push_back(-5);
		values.push_back(-10);
		values.push_back(20);
		values.push_back(7);

		StubMultiFilterOP sf( new StubMultiFilter );
		sf->set(values);
		protocols::filters::ReplicateFilter filter(sf,8);
		TS_ASSERT_EQUALS( filter.report_sm(*testpose_), (-2+-4+-1+-6+-5+-10+20+7)/8.0 );
		filter.median(true);
		TS_ASSERT_EQUALS( filter.report_sm(*testpose_), (-2+-4)/2.0 );

		protocols::filters::ReplicateFilter filter2(sf,8,/*upper*/ 1, /*lower*/ 3);
		TS_ASSERT_EQUALS( filter2.report_sm(*testpose_), (-2+-4+-1+7)/4.0 );
		filter2.median(true);
		TS_ASSERT_EQUALS( filter2.report_sm(*testpose_), (-2+-1)/2.0 );
	}

	void test_parsing() {
		basic::datacache::DataMap data;
		Filters_map filters;
		Movers_map movers;


		StubMultiFilterOP sf( new StubMultiFilter );
		utility::vector1<core::Real> values;
		values.push_back(-2);
		values.push_back(-4);
		values.push_back(-1);
		values.push_back(-6);
		values.push_back(-5);
		values.push_back(-10);
		values.push_back(20);
		values.push_back(7);
		values.push_back(-3);
		values.push_back(1);

		sf->set(values);

		filters["sf"] = sf;
		protocols::filters::ReplicateFilter  testfilter;
		TagCOP tag = tagptr_from_string("<ReplicateFilter name=test threshold=0 filter_name=sf replicates=10 upper_cut=0.2 lower_cut=4 />\n ");
		testfilter.parse_my_tag( tag, data, filters, movers, *testpose_ );

		TS_ASSERT_EQUALS( testfilter.report_sm( *testpose_), (-2+-1+-3+1)/4.0 );

		protocols::filters::ReplicateFilter  testfilter2;
		TagCOP tag2 = tagptr_from_string("<ReplicateFilter name=test2 median=1 threshold=0 filter_name=sf replicates=10 upper_cut=0.2 lower_cut=4 />\n ");
		testfilter2.parse_my_tag( tag2, data, filters, movers, *testpose_ );

		TS_ASSERT_EQUALS( testfilter2.report_sm( *testpose_), (-2+-1)/2.0 );
	}

};
