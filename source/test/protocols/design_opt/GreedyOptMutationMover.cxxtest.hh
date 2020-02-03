// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/minimization_packing/GreedyOptMutationMover.cxxtest.hh
/// @brief test for GreedyOptMutationMover; copy from GreenPacker (originally written to check GreenPacker bugs)
/// @author Steven Lewis smlewi@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>

// Project headers
#include <protocols/design_opt/GreedyOptMutationMover.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

// Package headers
#include <core/pose/Pose.hh>

#include <protocols/filters/BasicFilters.hh> //TrueFilter
#include <protocols/moves/NullMover.hh>

// Utility headers
#include <utility/vector1.hh>

// Numberic headers

// C++ headers


class GreedyOptMutationMoverTests : public CxxTest::TestSuite {

private:

public:

	void setUp() {
		protocols_init();
	}

	// @brief simple test: design on a small pose
	void test_GreedyOptMutationMover_ultrasimple() {

		core_init();
		core::pose::Pose pose(create_twores_1ubq_pose());
		//MQ is sequence of starting ubiquitin

		core::pack::task::TaskFactoryOP task_factory(new core::pack::task::TaskFactory());
		core::pack::task::operation::RestrictResidueToRepackingOP RRtROP( new core::pack::task::operation::RestrictResidueToRepacking());
		RRtROP->include_residue( 2 );
		task_factory->push_back( RRtROP );

		//trimming the design space to make the test MOAR FASTERER, because patience is for the weak
		core::pack::task::operation::RestrictAbsentCanonicalAASOP RACAAsOP( new core::pack::task::operation::RestrictAbsentCanonicalAAS());
		RACAAsOP->include_residue( 1 );
		utility::vector1<bool> const MOAR_FASTERER({1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
		RACAAsOP->keep_aas(MOAR_FASTERER);
		task_factory->push_back( RACAAsOP );

		protocols::design_opt::GreedyOptMutationMover GOMM;

		GOMM.task_factory( task_factory );
		GOMM.scorefxn(core::scoring::get_score_function());
		GOMM.skip_best_check(true); //for testing purposes, do not run the "check that all mutations are better than the WT was" check - we would have to have a test case where the WT was intentionally bad and I am not making that from a laptop on an airplane
		GOMM.add_filter(utility::pointer::make_shared<protocols::filters::TrueFilter>(), "high", 42);
		GOMM.relax_mover(utility::pointer::make_shared<protocols::moves::NullMover>()); //aint nobody got time for relax
		GOMM.apply( pose );

		//test that the sequence changed as expected
		TS_ASSERT_EQUALS("AQ", pose.sequence());

	}
}; //GreedyOptMutationMoverTests
