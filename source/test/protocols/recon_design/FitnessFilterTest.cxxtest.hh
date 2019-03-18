// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file test/protocols/recon_design/FitnessFilterTest.cxxtest.hh
/// @brief Unit test for FitnessFilter
/// @author Alex Sevy (alex.sevy@gmail.com)

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <protocols/recon_design/FitnessFilter.hh>

#include <basic/Tracer.hh>
#include <util/pose_funcs.hh>

static basic::Tracer TR("protocols.recon_design.FitnessFilterTest");

using namespace protocols::recon_design;
using namespace core::pose;
using namespace core::scoring::constraints;

class FitnessFilterTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		pose_ = create_trpcage_ideal_pose();
	}

	void tearDown() {
	}

	void test_FitnessFilter() {
		PoseOP pose1 = pose_.clone();
		PoseOP pose2 = pose_.clone();
		PoseOP pose3 = pose_.clone();
		utility::vector1<PoseOP> poses;
		poses.push_back( pose1 );
		poses.push_back( pose2 );
		poses.push_back( pose3 );

		FitnessFilter fitness;
		fitness.set_poses( poses );

		core::scoring::ScoreFunctionOP sfxn = fitness.sfxn();
		// each pose has fitness of 43.351558
		// total fitness is 130.054674
		fitness.threshold( 131 );
		TS_ASSERT( fitness.apply( *pose1 ) );
		fitness.threshold( 130 );
		TS_ASSERT( !fitness.apply( *pose1 ) );

		TS_ASSERT_DELTA( fitness.report_sm( *pose1 ), 130.054674, 0.0001 );

	}

private:
	core::pose::Pose pose_;
};
