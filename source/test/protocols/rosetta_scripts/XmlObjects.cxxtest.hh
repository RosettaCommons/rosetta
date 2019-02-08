// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/RosettaScriptsParser.cxxtest.hh
/// @brief Test suite for the RosettaScripts parser.
/// @details We were without a unit test for this until 6 May 2016!  This tests the xi:include functionality, but still isn't a comprehensive test.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <util/pose_funcs.hh>

#include <protocols/rosetta_scripts/XmlObjects.hh>

// Package headers
#include <core/pose/Pose.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/select/residue_selector/PhiSelector.hh>
#include <core/simple_metrics/test_classes.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/simple_filters/PoseCommentFilter.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

// Numberic headers

// C++ headers
#include <string>
#include <iostream>

static basic::Tracer TR("protocols.rosetta_scripts.XmlObjects.cxxtest");

using namespace protocols::rosetta_scripts;

////////////////////////////////////////////////////////////////////////
// Tests

class XmlObjectsTests : public CxxTest::TestSuite {


public:

	void setUp() {
		protocols_init();
		trpcage = create_trpcage_ideal_pose();
	}

	void test_string_instantiation() {

		TR << "Starting string instantiation test" << std::endl;

		std::string xml =
			"<SCOREFXNS>\n"
			"<ScoreFunction name=\"scorefxn\"/>\n"
			"</SCOREFXNS>\n"
			"<FILTERS>\n"
			"<PoseComment name=\"true\" />\n"
			"</FILTERS>\n"
			"<MOVERS>\n"
			"<MinMover name=\"min_mover\" bb=\"true\" chi=\"true\"/>\n"
			"</MOVERS>\n"
			"<TASKOPERATIONS>\n"
			"<PreventRepacking name=\"prevent_repacking\" />\n"
			"</TASKOPERATIONS>\n"
			"<RESIDUE_SELECTORS>\n"
			"<True name=\"true_sel\" />\n"
			"</RESIDUE_SELECTORS>\n"
			"<SIMPLE_METRICS>\n"
			"<TestRealMetric name=\"test_sm\" />\n"
			"</SIMPLE_METRICS>\n";

		XmlObjectsCOP objs = XmlObjects::create_from_string( xml );

		TR << "Getting ScoreFunction" << std::endl;
		try {
			core::scoring::ScoreFunctionOP scorefxn = objs->get_score_function("scorefxn");
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}

		TR << "Getting Filter" << std::endl;
		try {
			protocols::filters::FilterOP filter = objs->get_filter("true");
			protocols::simple_filters::PoseCommentOP casted = std::dynamic_pointer_cast< protocols::simple_filters::PoseComment >( filter );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}

		TR << "Getting Mover" << std::endl;
		try {
			protocols::moves::MoverOP mover = objs->get_mover("min_mover");
			protocols::minimization_packing::MinMoverOP casted = std::dynamic_pointer_cast< protocols::minimization_packing::MinMover >( mover );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}

		TR << "Getting TaskOperation" << std::endl;
		try {
			core::pack::task::operation::TaskOperationOP taskop = objs->get_task_operation("prevent_repacking");
			core::pack::task::operation::PreventRepackingOP casted = std::dynamic_pointer_cast< core::pack::task::operation::PreventRepacking >( taskop );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}

		TR << "Getting ResidueSelector" << std::endl;
		try {
			core::select::residue_selector::ResidueSelectorOP selector = objs->get_residue_selector("true_sel");
			core::select::residue_selector::TrueResidueSelectorOP casted = std::dynamic_pointer_cast< core::select::residue_selector::TrueResidueSelector >( selector );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}

		TR << "Getting SimpleMetric" << std::endl;
		try {
			core::simple_metrics::SimpleMetricOP metric = objs->get_simple_metric("test_sm");
			core::simple_metrics::TestRealMetricOP casted = std::dynamic_pointer_cast< core::simple_metrics::TestRealMetric >( metric );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}

		TR << "End string instantiation test" << std::endl;

	}


	// Originally this function had multiple semi-redundant tests
	// that just checked that it would work if you were missing
	// the name attribute. We need this test to be faster.

	void test_static_instantiation_score_function() {

		TR << "Starting static instantiation test" << std::endl;

		TR << "Getting ScoreFunction" << std::endl;
		try {
			core::scoring::ScoreFunctionOP scorefxn = XmlObjects::static_get_score_function(
				"<ScoreFunction name=\"scorefxn\"/>\n");
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TR << "End static instantiation score function test" << std::endl;
	}

	void dont_test_static_instantiation_residue_empty_sfxn() {

		try {
			core::scoring::ScoreFunctionOP scorefxn = XmlObjects::static_get_score_function(
				"<ScoreFunction />\n");
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TR << "End static instantiation empty sfxn test" << std::endl;
	}

	void test_static_instantiation_filter() {


		TR << "Getting Filter" << std::endl;
		try {
			protocols::filters::FilterOP filter = XmlObjects::static_get_filter(
				"<PoseComment name=\"true\" />\n");
			protocols::simple_filters::PoseCommentOP casted = std::dynamic_pointer_cast< protocols::simple_filters::PoseComment >( filter );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TR << "End static instantiation filter test" << std::endl;
	}

	void dont_test_static_instantiation_empty_filter() {

		try {
			protocols::filters::FilterOP filter = XmlObjects::static_get_filter(
				"<PoseComment />\n");
			protocols::simple_filters::PoseCommentOP casted = std::dynamic_pointer_cast< protocols::simple_filters::PoseComment >( filter );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TR << "End static instantiation empty filter test" << std::endl;
	}

	void test_static_instantiation_mover() {

		TR << "Getting Mover" << std::endl;
		try {
			protocols::moves::MoverOP mover = XmlObjects::static_get_mover(
				"<MinMover name=\"min_mover\" bb=\"true\" chi=\"true\"/>\n");
			protocols::minimization_packing::MinMoverOP casted = std::dynamic_pointer_cast< protocols::minimization_packing::MinMover >( mover );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TR << "End static instantiation mover test" << std::endl;
	}

	void dont_test_static_instantiation_mover_noname() {

		try {
			protocols::moves::MoverOP mover = XmlObjects::static_get_mover(
				"<MinMover bb=\"true\" chi=\"true\"/>\n");
			protocols::minimization_packing::MinMoverOP casted = std::dynamic_pointer_cast< protocols::minimization_packing::MinMover >( mover );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TR << "End static instantiation mover noname test" << std::endl;
	}

	void test_static_instantiation_taskop() {


		TR << "Getting TaskOperation" << std::endl;
		try {
			core::pack::task::operation::TaskOperationOP taskop = XmlObjects::static_get_task_operation(
				"<PreventRepacking name=\"prevent_repacking\" />\n");
			core::pack::task::operation::PreventRepackingOP casted = std::dynamic_pointer_cast< core::pack::task::operation::PreventRepacking >( taskop );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TR << "End static instantiation taskop test" << std::endl;
	}

	void dont_test_static_instantiation_taskop_noname() {

		try {
			core::pack::task::operation::TaskOperationOP taskop = XmlObjects::static_get_task_operation(
				"<PreventRepacking />\n");
			core::pack::task::operation::PreventRepackingOP casted = std::dynamic_pointer_cast< core::pack::task::operation::PreventRepacking >( taskop );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TR << "End static instantiation taskop noname test" << std::endl;
	}

	void test_static_instantiation_residue_selector() {


		TR << "Getting ResidueSelector" << std::endl;
		try {
			core::select::residue_selector::ResidueSelectorOP selector = XmlObjects::static_get_residue_selector(
				"<True name=\"true_sel\" />\n");
			core::select::residue_selector::TrueResidueSelectorOP casted = std::dynamic_pointer_cast< core::select::residue_selector::TrueResidueSelector >( selector );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TR << "End static instantiation residue selector test" << std::endl;
	}

	void dont_test_static_instantiation_residue_selector_noname() {

		try {
			core::select::residue_selector::ResidueSelectorOP selector = XmlObjects::static_get_residue_selector(
				"<True />\n");
			core::select::residue_selector::TrueResidueSelectorOP casted = std::dynamic_pointer_cast< core::select::residue_selector::TrueResidueSelector >( selector );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}


		TR << "End static instantiation residue selector noname test" << std::endl;
	}

	void test_static_instantiation_simple_metric() {


		TR << "Getting SimpleMetric" << std::endl;
		try {
			core::simple_metrics::SimpleMetricOP metric = XmlObjects::static_get_simple_metric(
				"<TestRealMetric name=\"test_sm\" />\n");
			core::simple_metrics::TestRealMetricOP casted = std::dynamic_pointer_cast< core::simple_metrics::TestRealMetric >( metric );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TR << "End static instantiation residue selector test" << std::endl;
	}

	void dont_test_static_instantiation_simple_metric_noname() {

		try {
			core::simple_metrics::SimpleMetricOP metric = XmlObjects::static_get_simple_metric(
				"<TestRealMetric />\n");
			core::simple_metrics::TestRealMetricOP casted = std::dynamic_pointer_cast< core::simple_metrics::TestRealMetric >( metric );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}


		TR << "End static instantiation residue selector noname test" << std::endl;
	}

	void test_file_instantiation() {

		TR << "Starting file instantiation test" << std::endl;

		XmlObjectsCOP objs = XmlObjects::create_from_file( "protocols/rosetta_scripts/modern_test1.xml" );

		try {
			core::select::residue_selector::ResidueSelectorOP selector = objs->get_residue_selector( "posPhi" );
			core::select::residue_selector::PhiSelectorOP casted = std::dynamic_pointer_cast< core::select::residue_selector::PhiSelector >( selector );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}


		TR << "End file instantiation test" << std::endl;
	}


	void test_includes() {

		TR << "Starting includes test" << std::endl;

		XmlObjectsCOP objs = XmlObjects::create_from_string( "<xi:include href=\"protocols/rosetta_scripts/modern_include2.xml\" / >" );

		try {
			core::select::residue_selector::ResidueSelectorOP selector = objs->get_residue_selector( "posPhi" );
			core::select::residue_selector::PhiSelectorOP casted = std::dynamic_pointer_cast< core::select::residue_selector::PhiSelector >( selector );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}

		TR << "End includes test" << std::endl;

	}




	void test_parsed_protocol() {

		TR << "Starting parsed protocol test" << std::endl;


		XmlObjectsCOP objs = XmlObjects::create_from_file( "protocols/rosetta_scripts/modern_test1.xml" );
		try {
			protocols::moves::MoverOP mover = objs->get_mover( "ParsedProtocol" );
			protocols::rosetta_scripts::ParsedProtocolOP casted = std::dynamic_pointer_cast< protocols::rosetta_scripts::ParsedProtocol >( mover );
			(void)casted;
		} catch (utility::excn::Exception const & e ) {
			TR << e.msg() << std::endl;
			TS_ASSERT( false );
		}


		TR << "End parsed protocol test" << std::endl;
	}


private:
	core::pose::Pose trpcage;


};
