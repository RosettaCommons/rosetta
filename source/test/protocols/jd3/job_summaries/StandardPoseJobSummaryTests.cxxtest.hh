// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/jd3/job_summaries/StandardPoseJobSummaryTests.cxxtest.hh
/// @brief  Tests for the StandardPoseJobSummary
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/jd3/job_summaries/StandardPoseJobSummary.hh>
#include <core/simple_metrics/test_classes.hh>
#include <core/simple_metrics/SimpleMetricData.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("StandardPoseJobSummaryTests");

using namespace core::scoring;
using namespace core::simple_metrics;
using namespace protocols::jd3::job_summaries;

class StandardPoseJobSummaryTests : public CxxTest::TestSuite {
	//Define Variables

public:

	core::pose::Pose pose;

	void setUp(){
		core_init();
		core::import_pose::pose_from_file(pose, "protocols/antibody/aho_with_antigen.pdb", core::import_pose::PDB_file);
	}

	void tearDown(){

	}

	void test_basic_pose_job_summary(){
		ScoreFunctionOP scorefxn = get_score_function();
		core::Real pose_e = scorefxn->score(pose);

		TestRealMetricOP real_metric = TestRealMetricOP( new TestRealMetric());
		TestStringMetricOP str_metric = TestStringMetricOP( new TestStringMetric());

		StandardPoseJobSummaryOP job_summary = StandardPoseJobSummaryOP( new StandardPoseJobSummary());
		job_summary->set_energy(10);

		TS_ASSERT_EQUALS(job_summary->energy(), 10);

		real_metric->apply(pose);
		str_metric->apply(pose);

		job_summary->extract_summary(pose);

		TS_ASSERT_DELTA(job_summary->energy(), pose_e, .01)

			SimpleMetricDataCOP metric_data = job_summary->metric_data();
		core::Real real_metric_value;
		std::string str_metric_value;

		metric_data->get_value("SomeReal", real_metric_value);
		metric_data->get_value("SomeString", str_metric_value);

		TS_ASSERT_DELTA(real_metric_value, 1.0, .001);
		TS_ASSERT_EQUALS(str_metric_value, "TESTING");
	}






};
