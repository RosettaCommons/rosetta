// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/filters/BasicFilters.cxxtest.hh
/// @brief  test for basic filters
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/rosettascripts.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>

// Utility Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.filters.BasicFilters.cxxtest.hh");

// --------------- Test Class --------------- //

class StochasticFilterTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP testpose_;
public:

	void setUp() {
		core_init();

		testpose_ = create_twores_1ubq_poseop(); // Identity doesn't matter
	}

	void tearDown() {
	}

	// This is a test of the direct application of StochasticFilter
	void test_apply() {
		utility::vector1< protocols::filters::StochasticFilter > filters = { {1.0}, {0.75}, {0.25}, {0.0} };
		utility::vector1< core::Size > mins = { 100, 67, 17, 0 };
		utility::vector1< core::Size > maxes = { 100, 83, 33, 0 };
		utility::vector1< core::Size > counts = { 0, 0, 0, 0 };

		for ( core::Size ii(1); ii <= 100; ++ii ) { // 100 replicates
			for ( core::Size jj(1); jj <= filters.size(); ++jj ) {
				counts[jj] += filters[jj].apply(*testpose_);
			}
		}

		for ( core::Size jj(1); jj <= filters.size(); ++jj ) {
			TR << "For StochasticFilter with confidence " << filters[jj].confidence() << " got " << counts[jj] << " true values. Expected between " << mins[jj] << " and " << maxes[jj] << std::endl;
			TS_ASSERT( counts[jj] >= mins[ jj ] );
			TS_ASSERT( counts[jj] <= maxes[ jj ] );
		}

	}

	void test_subfilter_apply_true() {
		using namespace protocols::filters;

		utility::vector1< core::Real > confidences = { 1.0, 0.75, 0.25, 0.0 };
		utility::vector1< core::Size > num_trues;
		utility::vector1< core::Size > counts;
		utility::vector1< core::Size > report_sm_counts;

		for ( core::Real confidence: confidences ) {
			StubFilterOP subfilter( new StubFilter(false, 3.5) );

			FilterOP stochastic_filter( new StochasticFilter( confidence, subfilter, /*run_subfilter_on*/ true ) );

			core::Size num_true=0;
			for ( core::Size ii(1); ii <= 100; ++ii ) { // 100 replicates
				num_true += stochastic_filter->apply(*testpose_);
				TS_ASSERT_EQUALS( stochastic_filter->report_sm(*testpose_), 3.5 ); // Should always call the subfilter
			}
			num_trues.push_back( num_true );
			counts.push_back( subfilter->num_apply );
			report_sm_counts.push_back( subfilter->num_report_sm );
		}

		for ( core::Size jj(1); jj <= confidences.size(); ++jj ) {
			TR << "For a confidence of " << confidences[jj] << " got " << num_trues[jj] << " true results." << std::endl;
			TS_ASSERT_EQUALS( num_trues[jj], 0 ); // Run on true, but returns false

			TR << "For a confidence of " << confidences[jj] << " got " << counts[jj] << " sub-filter applications." << std::endl;
			TS_ASSERT( counts[jj] >= (confidences[jj]-0.08)*100 );
			TS_ASSERT( counts[jj] <= (confidences[jj]+0.08)*100 );

			TR << "For a confidence of " << confidences[jj] << " got " << report_sm_counts[jj] << " sub-filter report_sm." << std::endl;
			TS_ASSERT_EQUALS( report_sm_counts[jj], 100 );
		}
	}

	void test_subfilter_apply_false() {
		using namespace protocols::filters;

		utility::vector1< core::Real > confidences = { 1.0, 0.75, 0.25, 0.0 };
		utility::vector1< core::Size > num_trues;
		utility::vector1< core::Size > counts;
		utility::vector1< core::Size > report_sm_counts;

		for ( core::Real confidence: confidences ) {
			StubFilterOP subfilter( new StubFilter(true, 3.5) );

			FilterOP stochastic_filter( new StochasticFilter( confidence, subfilter, /*run_subfilter_on*/ false ) );

			core::Size num_true=0;
			for ( core::Size ii(1); ii <= 100; ++ii ) { // 100 replicates
				num_true += stochastic_filter->apply(*testpose_);
				TS_ASSERT_EQUALS( stochastic_filter->report_sm(*testpose_), 3.5 ); // Should always call the subfilter
			}
			num_trues.push_back( num_true );
			counts.push_back( subfilter->num_apply );
			report_sm_counts.push_back( subfilter->num_report_sm );
		}

		for ( core::Size jj(1); jj <= confidences.size(); ++jj ) {
			TR << "For a confidence of " << confidences[jj] << " got " << num_trues[jj] << " true results. (run on false)" << std::endl;
			TS_ASSERT_EQUALS( num_trues[jj], 100 ); // Run on false, but then returns true

			TR << "For a confidence of " << confidences[jj] << " got " << counts[jj] << " sub-filter applications. (run on false)" << std::endl;
			TS_ASSERT( counts[jj] <= (1.0-(confidences[jj]-0.08))*100 );
			TS_ASSERT( counts[jj] >= (1.0-(confidences[jj]+0.08))*100 );

			TR << "For a confidence of " << confidences[jj] << " got " << report_sm_counts[jj] << " sub-filter report_sm. (run on false)" << std::endl;
			TS_ASSERT_EQUALS( report_sm_counts[jj], 100 );
		}
	}
};


class IfThenFilterTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP testpose_;
public:

	void setUp() {
		core_init();

		testpose_ = create_twores_1ubq_poseop(); // Identity doesn't matter
	}

	void tearDown() {
	}

	void test_elseonly() {
		protocols::filters::IfThenFilter filter;
		filter.threshold( 0 );
		filter.set_else( 0 , -10 );
		TS_ASSERT( filter.apply(*testpose_) );
		TS_ASSERT_EQUALS( filter.report_sm(*testpose_), -10 );

		filter.set_else( 0 , 11.25 );
		TS_ASSERT( ! filter.apply(*testpose_) );
		TS_ASSERT_EQUALS( filter.report_sm(*testpose_), 11.25 );

		protocols::filters::FilterOP sf5( new StubFilter( true, -5 ) );
		filter.set_else( sf5 , 11.25 );
		TS_ASSERT( filter.apply(*testpose_) );
		TS_ASSERT_EQUALS( filter.report_sm(*testpose_), -5 );

		protocols::filters::FilterOP sf11( new StubFilter( false, 10.875) );
		filter.set_else( sf11 , 11.25 );
		TS_ASSERT( ! filter.apply(*testpose_) );
		TS_ASSERT_EQUALS( filter.report_sm(*testpose_), 10.875 );

		protocols::filters::FilterOP sf12( new StubFilter( true, 11.125) );
		filter.set_else( sf12 , 11.25 );
		TS_ASSERT( ! filter.apply(*testpose_) );
		TS_ASSERT_EQUALS( filter.report_sm(*testpose_), 11.125 );
	}

	void test_oneif() {
		protocols::filters::FilterOP sfT( new StubFilter( true, 222) );
		protocols::filters::FilterOP sfF( new StubFilter( false, 225) );
		protocols::filters::FilterOP sfV( new StubFilter( true, 335) );

		protocols::filters::IfThenFilter filter1;
		filter1.threshold( 0 );
		filter1.set_else( 0 , -111 );
		filter1.add_condition( sfT, 0, 333 );
		TS_ASSERT_EQUALS( filter1.report_sm(*testpose_), 333 );

		protocols::filters::IfThenFilter filter2;
		filter2.threshold( 0 );
		filter2.set_else( 0 , -111 );
		filter2.add_condition( sfF, 0, 333 );
		TS_ASSERT_EQUALS( filter2.report_sm(*testpose_), -111 );

		protocols::filters::IfThenFilter filter3;
		filter3.threshold( 0 );
		filter3.set_else( 0 , -111 );
		filter3.add_condition( sfT, sfV, 333 );
		TS_ASSERT_EQUALS( filter3.report_sm(*testpose_), 335 );

		protocols::filters::IfThenFilter filter4;
		filter4.threshold( 0 );
		filter4.set_else( 0 , -111 );
		filter4.add_condition( sfF, sfV, 333 );
		TS_ASSERT_EQUALS( filter4.report_sm(*testpose_), -111 );
	}

	void test_multi() {
		protocols::filters::FilterOP sfF1( new StubFilter( false, 111) );
		protocols::filters::FilterOP sfF2( new StubFilter( false, 222) );
		protocols::filters::FilterOP sfF3( new StubFilter( false, 333) );
		protocols::filters::FilterOP sfT( new StubFilter( true, 225) );

		protocols::filters::IfThenFilter filter;
		filter.threshold( 0 );
		filter.set_else( 0 , -10 );

		filter.add_condition( sfF1, 0, 1110 );
		filter.add_condition( sfF2, 0, 2220 );
		filter.add_condition( sfT, 0, 9999 );
		filter.add_condition( sfF3, 0, 3330 );

		TS_ASSERT_EQUALS( filter.report_sm(*testpose_), 9999 );

		protocols::filters::IfThenFilter filter2;
		filter2.threshold( 0 );
		filter2.set_else( 0 , -10 );

		filter2.add_condition( sfF1, 0, 11110 );
		filter2.add_condition( sfF2, 0, 2220 );
		filter2.add_condition( sfF3, 0, 3330 );

		TS_ASSERT_EQUALS( filter2.report_sm(*testpose_), -10 );
	}

	void test_parsing() {
		basic::datacache::DataMap data;
		Filters_map filters;
		Movers_map movers;

		filters["sfF1"] = protocols::filters::FilterOP( new StubFilter( false, 1 ) );
		filters["sfF2"] = protocols::filters::FilterOP( new StubFilter( false, 2 ) );
		filters["sfF3"] = protocols::filters::FilterOP( new StubFilter( false, 3 ) );
		filters["sfT10"] = protocols::filters::FilterOP( new StubFilter( true, 10 ) );
		filters["sfT20"] = protocols::filters::FilterOP( new StubFilter( true, 20 ) );
		filters["sfT99"] = protocols::filters::FilterOP( new StubFilter( true, 99 ) );

		protocols::filters::IfThenFilter  testfilter;
		TagCOP tag = tagptr_from_string("<IfThenFilter name=test threshold=0 >\n"
			"    <IF testfilter=sfF1 valuefilter=sfT10 />\n"
			"    <ELIF testfilter=sfF2 value=900 />\n"
			"    <ELIF testfilter=sfT20 valuefilter=sfF3 />\n"
			"    <ELSE valuefilter=sfT99 />\n"
			" </IfThenFilter>\n");

		testfilter.parse_my_tag( tag, data, filters, movers, *testpose_ );

		TS_ASSERT_EQUALS( testfilter.report_sm( *testpose_), 3 );
	}

	void test_weights() {
		protocols::filters::IfThenFilter filter;
		filter.threshold( 0 );
		filter.set_else( 0 , -10, 2 );
		TS_ASSERT_EQUALS( filter.report_sm(*testpose_), -20 );

		protocols::filters::FilterOP sf5( new StubFilter( true, -5 ) );
		filter.set_else( sf5 , 11.25, 3 );
		TS_ASSERT_EQUALS( filter.report_sm(*testpose_), -15 );

		protocols::filters::FilterOP sfT( new StubFilter( true, 222) );
		protocols::filters::FilterOP sfF( new StubFilter( false, 225) );
		protocols::filters::FilterOP sfV( new StubFilter( true, 335) );

		protocols::filters::IfThenFilter filter1;
		filter1.threshold( 0 );
		filter1.set_else( 0 , -111 );
		filter1.add_condition( sfT, 0, 333, false, 3 );
		TS_ASSERT_EQUALS( filter1.report_sm(*testpose_), 999 );

		protocols::filters::IfThenFilter filter2;
		filter2.threshold( 0 );
		filter2.set_else( 0 , -112 );
		filter2.add_condition( sfT, sfV, 333, false, 2 );
		TS_ASSERT_EQUALS( filter2.report_sm(*testpose_), 670 );
	}

	void test_floor() {
		protocols::filters::IfThenFilter filter;
		filter.threshold( -5 );
		filter.set_else( 0 , -10, 2 );
		TS_ASSERT( filter.apply(*testpose_) );

		filter.set_lower_threshold( true );
		TS_ASSERT( ! filter.apply(*testpose_) );

		filter.threshold( -30 );
		TS_ASSERT( filter.apply(*testpose_) );
	}

	void test_invert() {
		protocols::filters::FilterOP sfF1( new StubFilter( false, 111) );
		protocols::filters::FilterOP sfF2( new StubFilter( false, 222) );
		protocols::filters::FilterOP sfF3( new StubFilter( false, 333) );
		protocols::filters::FilterOP sfT( new StubFilter( true, 225) );

		protocols::filters::IfThenFilter filter;
		filter.threshold( 0 );
		filter.set_else( 0 , -10 );

		filter.add_condition( sfF1, 0, 1110, false );
		filter.add_condition( sfT, 0, 2220, true );
		filter.add_condition( sfF2, 0, 9999, true );
		filter.add_condition( sfF3, 0, 3330 );

		TS_ASSERT_EQUALS( filter.report_sm(*testpose_), 9999 );
	}

	void test_parse_invertweights() {
		basic::datacache::DataMap data;
		Filters_map filters;
		Movers_map movers;

		filters["sfF1"] = protocols::filters::FilterOP( new StubFilter( false, 1 ) );
		filters["sfF2"] = protocols::filters::FilterOP( new StubFilter( false, 2 ) );
		filters["sfF3"] = protocols::filters::FilterOP( new StubFilter( false, 3 ) );
		filters["sfT10"] = protocols::filters::FilterOP( new StubFilter( true, 10 ) );
		filters["sfT20"] = protocols::filters::FilterOP( new StubFilter( true, 20 ) );
		filters["sfT99"] = protocols::filters::FilterOP( new StubFilter( true, 99 ) );

		protocols::filters::IfThenFilter  testfilter;
		TagCOP tag = tagptr_from_string("<IfThenFilter name=test threshold=2 lower_threshold=1>\n"
			"    <IF testfilter=sfF1 valuefilter=sfT10 weight=60/>\n"
			"    <ELIF testfilter=sfF2 inverttest=1 valuefilter=sfF3 weight=2/>\n"
			"    <ELSE valuefilter=sfT99 weight=5/>\n"
			" </IfThenFilter>\n");

		testfilter.parse_my_tag( tag, data, filters, movers, *testpose_ );

		TS_ASSERT_EQUALS( testfilter.report_sm( *testpose_), 6 );
		TS_ASSERT( testfilter.apply(*testpose_) );
	}

};
