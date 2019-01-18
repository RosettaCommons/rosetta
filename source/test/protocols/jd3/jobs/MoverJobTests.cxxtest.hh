// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/jd3/jobs/MoverJobTests.cxxtest.hh
/// @brief  Tests for the MoverJob
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
#include <core/simple_metrics/test_classes.hh>
#include <core/simple_metrics/SimpleMetricData.hh>
#include <core/conformation/Conformation.hh>

#include <protocols/carbohydrates/CreateGlycanSequonMover.hh>
#include <protocols/carbohydrates/SimpleGlycosylateMover.hh>
#include <protocols/jd3/job_summaries/StandardPoseJobSummary.hh>
#include <protocols/jd3/jobs/MoverJob.hh>
#include <protocols/jd3/job_results/PoseJobResult.hh>
#include <protocols/moves/MoverContainer.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("MoverJobTests");

using namespace core::scoring;
using namespace core::simple_metrics;
using namespace protocols::jd3::job_summaries;
using namespace protocols::jd3::job_results;
using namespace protocols::jd3::jobs;
using namespace protocols::carbohydrates;
using namespace protocols::moves;

class MoverJobTests : public CxxTest::TestSuite {
	//Define Variables

public:
	void setUp(){
		core_init_with_additional_options("-include_sugars");
	}

	void tearDown(){

	}

	void test_basic_mover_job(){
		core::pose::PoseOP pose = core::import_pose::pose_from_file( "protocols/antibody/aho_with_antigen.pdb", false, core::import_pose::PDB_file);
		ScoreFunctionOP scorefxn = get_score_function();
		//core::Real pose_e = scorefxn->score(*pose);

		TestRealMetricOP real_metric = TestRealMetricOP( new TestRealMetric());
		TestStringMetricOP str_metric = TestStringMetricOP( new TestStringMetric());
		utility::vector1< SimpleMetricCOP > metrics;
		metrics.push_back(real_metric);
		metrics.push_back(str_metric);

		//55
		SimpleGlycosylateMoverOP glycosylate = SimpleGlycosylateMoverOP( new SimpleGlycosylateMover());
		glycosylate->set_position(55);
		glycosylate->set_glycosylation("man5");

		CreateGlycanSequonMoverOP sequon_creator = CreateGlycanSequonMoverOP( new CreateGlycanSequonMover());
		sequon_creator->set_pack_neighbors(false);
		sequon_creator->set_glycosylation_position(55, *pose);

		MoverJobOP mover_job = MoverJobOP( new MoverJob());
		SequenceMoverOP seq_mover = SequenceMoverOP( new SequenceMover(sequon_creator, glycosylate));

		mover_job->set_mover(seq_mover);
		mover_job->add_metric(real_metric, "job1_set1_");
		mover_job->add_metric(str_metric, "job1_set1_");

		utility::vector1< core::simple_metrics::SimpleMetricCOP > set2_metrics;
		set2_metrics.push_back(real_metric);
		set2_metrics.push_back(str_metric);

		mover_job->add_metrics(set2_metrics, "job1_set2_");

		mover_job->pose(pose);

		protocols::jd3::CompletedJobOutput output = mover_job->run();
		StandardPoseJobSummaryOP job_summary = utility::pointer::dynamic_pointer_cast< StandardPoseJobSummary > ( output.job_results[1].first );
		PoseJobResultOP pose_job_result = utility::pointer::dynamic_pointer_cast< PoseJobResult > ( output.job_results[1].second );

		TS_ASSERT(pose_job_result->pose()->conformation().contains_carbohydrate_residues());

		SimpleMetricDataCOP metric_data = job_summary->metric_data();
		core::Real real_metric_value;
		std::string str_metric_value;

		metric_data->get_value("job1_set1_SomeReal", real_metric_value);
		metric_data->get_value("job1_set1_SomeString", str_metric_value);

		TS_ASSERT_DELTA(real_metric_value, 1.0, .001);
		TS_ASSERT_EQUALS(str_metric_value, "TESTING");

		metric_data->get_value("job1_set2_SomeReal", real_metric_value);
		metric_data->get_value("job1_set2_SomeString", str_metric_value);

		TS_ASSERT_DELTA(real_metric_value, 1.0, .001);
		TS_ASSERT_EQUALS(str_metric_value, "TESTING");

	}

	void test_seq_mover_job(){
		core::pose::PoseOP pose = core::import_pose::pose_from_file( "protocols/antibody/aho_with_antigen.pdb", false, core::import_pose::PDB_file);

		ScoreFunctionOP scorefxn = get_score_function();
		//core::Real pose_e = scorefxn->score(*pose);

		TS_ASSERT(! pose->conformation().contains_carbohydrate_residues());

		TestRealMetricOP real_metric = TestRealMetricOP( new TestRealMetric());
		TestStringMetricOP str_metric = TestStringMetricOP( new TestStringMetric());
		utility::vector1< SimpleMetricCOP > metrics;
		metrics.push_back(real_metric);
		metrics.push_back(str_metric);

		//55
		SimpleGlycosylateMoverOP glycosylate = SimpleGlycosylateMoverOP( new SimpleGlycosylateMover());
		glycosylate->set_position(55);
		glycosylate->set_glycosylation("man5");

		CreateGlycanSequonMoverOP sequon_creator = CreateGlycanSequonMoverOP( new CreateGlycanSequonMover());
		sequon_creator->set_pack_neighbors(false);
		sequon_creator->set_glycosylation_position(55, *pose);

		MoverJobOP mover_job = MoverJobOP( new MoverJob());

		mover_job->add_mover(sequon_creator);
		mover_job->add_mover(glycosylate);
		mover_job->add_metric(real_metric, "job1_set1_");
		mover_job->add_metric(str_metric, "job1_set1_");

		utility::vector1< core::simple_metrics::SimpleMetricCOP > set2_metrics;
		set2_metrics.push_back(real_metric);
		set2_metrics.push_back(str_metric);

		mover_job->add_metrics(set2_metrics, "job1_set2_");

		mover_job->pose(pose);

		protocols::jd3::CompletedJobOutput output = mover_job->run();
		StandardPoseJobSummaryOP job_summary = utility::pointer::dynamic_pointer_cast< StandardPoseJobSummary > ( output.job_results[1].first );
		PoseJobResultOP pose_job_result = utility::pointer::dynamic_pointer_cast< PoseJobResult > ( output.job_results[1].second );

		TS_ASSERT(pose_job_result->pose()->conformation().contains_carbohydrate_residues());

		SimpleMetricDataCOP metric_data = job_summary->metric_data();
		core::Real real_metric_value;
		std::string str_metric_value;

		metric_data->get_value("job1_set1_SomeReal", real_metric_value);
		metric_data->get_value("job1_set1_SomeString", str_metric_value);

		TS_ASSERT_DELTA(real_metric_value, 1.0, .001);
		TS_ASSERT_EQUALS(str_metric_value, "TESTING");

		metric_data->get_value("job1_set2_SomeReal", real_metric_value);
		metric_data->get_value("job1_set2_SomeString", str_metric_value);

		TS_ASSERT_DELTA(real_metric_value, 1.0, .001);
		TS_ASSERT_EQUALS(str_metric_value, "TESTING");

	}

};
