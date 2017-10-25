// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/filters/CalculatorFilter.cxxtest.hh
/// @brief  test for the Calculator filter
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/rosettascripts.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/filters/CalculatorFilter.hh>

// Utility Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.filters.CalculatorFilter.cxxtest.hh");

// --------------- Test Class --------------- //

class CalculatorFilterTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP testpose_;
public:

	void setUp() {
		core_init();

		testpose_ = create_twores_1ubq_poseop(); // Identity doesn't matter
	}

	void tearDown() {
	}

	void test_calculatorfilter() {
		protocols::filters::CalculatorFilter cf(" t1 = exp(-E1/kT); t2 = exp(-E2/kT); t1/( t1 + t2 ) " );

		cf.add_constant("kT", 0.6 );
		cf.add_filter("E1", protocols::filters::FilterOP( new StubFilter( true, -2) ) );
		cf.add_filter("E2", protocols::filters::FilterOP( new StubFilter( true, -1) ) );

		//default 0 threshold
		TS_ASSERT( ! cf.apply(*testpose_) );
		TS_ASSERT_DELTA( cf.report_sm(*testpose_), 0.84113089511, 0.0001  );
		cf.threshold( 2.0 );
		TS_ASSERT( cf.apply(*testpose_) );
	}

	void test_parsing() {
		basic::datacache::DataMap data;
		Filters_map filters;
		Movers_map movers;

		StubFilterOP sf1( new StubFilter( true, -1) );
		StubFilterOP sf2( new StubFilter( true, -2) );
		StubFilterOP sf3( new StubFilter( true, -3) );

		filters["alpha"] = sf1;
		filters["beta"] = sf2;
		filters["delta"] = sf3;

		protocols::filters::CalculatorFilter  testfilter;
		TagCOP tag = tagptr_from_string("<CalculatorFilter name=test threshold=0 equation=\"(min(a, b, c)/(a*b*c-d)) + e\" >\n "
			"<var name=a filter=alpha />\n"
			"<Var name=b filter_name=beta />\n"
			"<Var name=c filter=delta />\n"
			"<VAR name=d value=-4.0 />\n"
			"<VAR name=e reported=ten />\n"
			"</CalculatorFilter>\n" );

		testfilter.parse_my_tag( tag, data, filters, movers, *testpose_ );
		setPoseExtraScore(*testpose_, "ten", 10.0);


		TS_ASSERT_EQUALS( testfilter.report_sm( *testpose_), 11.5 );
	}

};
