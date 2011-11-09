// -*- mode:c++;tab-width:2;indent-tabs-mode:nil;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rotamer_recovery/RotamerRecovery.cxxtest.hh
/// @brief  Test RotamerRecovery class
/// @author Matthew O'Meara (mattjomeara@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <util/pose_funcs.hh>

// Unit Headers
#include <protocols/rotamer_recovery/RRProtocol.hh>
#include <protocols/rotamer_recovery/RRProtocolRTMin.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRComparerAutomorphicRMSD.hh>
#include <protocols/rotamer_recovery/RotamerRecovery.hh>
#include <protocols/rotamer_recovery/RotamerRecoveryFactory.hh>

// Project Headers
#include <test/core/init_util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// C++ Headers
#include <iostream>

//Auto Headers
#include <utility/vector1.hh>



static basic::Tracer TR("protocols.rotamer_recovery.RotamerRecovery.cxxtest");

class RotamerRecoveryTests : public CxxTest::TestSuite {

public:

	void
	setUp() {
		using core::scoring::getScoreFunction;
		using core::pack::task::TaskFactory;
		using core::pack::task::operation::RestrictToRepacking;

		core_init();
		pose_1ten_ = fullatom_pose_from_string( pdb_string_1ten() );
		score_function_ = getScoreFunction();

		TaskFactory task_factory;

		task_factory.push_back( new RestrictToRepacking );
		packer_task_1ten_ = task_factory.create_task_and_apply_taskoperations( pose_1ten_ );

	}

	void test_RotamerRecovery_main() {
		do_test_RotamerRecovery_automatic_construction();
	}

	void
	do_test_RotamerRecovery_automatic_construction() {

		using std::endl;
		using core::Real;
		using protocols::rotamer_recovery::RotamerRecoveryOP;
    using protocols::rotamer_recovery::RotamerRecoveryFactory;
		score_function_->setup_for_scoring(pose_1ten_);

    RotamerRecoveryFactory* factory(RotamerRecoveryFactory::get_instance());
		{
      RotamerRecoveryOP rr(
        factory->get_rotamer_recovery(
          "RRProtocolRTMin", "RRComparerRotBins", "RRReporterSimple"));

			rr->run( pose_1ten_,*score_function_,*packer_task_1ten_);
			rr->show( TR );
			TS_ASSERT_DELTA( rr->recovery_rate(), Real(64)/Real(89) , .001 );
		}

		{
      RotamerRecoveryOP rr(
        factory->get_rotamer_recovery(
          "RRProtocolRTMin", "RRComparerAutomorphicRMSD", "RRReporterSimple"));
			rr->run( pose_1ten_,*score_function_,*packer_task_1ten_);
			rr->show( TR );
			TS_ASSERT_DELTA( rr->recovery_rate(), Real(22)/Real(89) , .001 );
		}


	}

private:
	core::pose::Pose pose_1ten_;
	core::scoring::ScoreFunctionOP score_function_;
	core::pack::task::PackerTaskOP packer_task_1ten_;

};
