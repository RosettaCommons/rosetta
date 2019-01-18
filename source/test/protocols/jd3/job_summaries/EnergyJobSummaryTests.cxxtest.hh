// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/jd3/job_summaries/EnergyJobSummaryTests.cxxtest.hh
/// @brief  Tests for the EnergyJobSummary.
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
#include <protocols/jd3/job_summaries/EnergyJobSummary.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("EnergyJobSummaryTests");

using namespace core::scoring;
using namespace protocols::jd3::job_summaries;

class EnergyJobSummaryTests : public CxxTest::TestSuite {
	//Define Variables

public:
	core::pose::Pose pose;

	void setUp(){
		core_init();
		core::import_pose::pose_from_file(pose, "protocols/antibody/aho_with_antigen.pdb", core::import_pose::PDB_file);
	}

	void tearDown(){

	}

	void test_basic_energy_job_summary() {
		ScoreFunctionOP scorefxn = get_score_function();
		core::Real pose_e = scorefxn->score(pose);

		EnergyJobSummaryOP job_summary = EnergyJobSummaryOP( new EnergyJobSummary(10));
		TS_ASSERT_EQUALS(job_summary->energy(), 10);

		job_summary->extract_energy(pose);
		TS_ASSERT_DELTA(job_summary->energy(), pose_e, .01);

		job_summary->energy(10);
		TS_ASSERT_EQUALS(job_summary->energy(), 10);

	}






};
