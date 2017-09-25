// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/OperateOnResidueSubset.cxxtest.hh
/// @brief  test suite for core::select::residue_selector::OperateOnResidueSubset
/// @author Steven Lewis smlewi@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>

// Package headers
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/FillAUTOTaskOperation.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::pack::task;
using namespace core::pack::task::operation;

class FillAUTOTaskOperationTests : public CxxTest::TestSuite {

public:

	core::pose::PoseOP poseop;

	void setUp() {
		core_init();

		if ( !poseop ) poseop = create_twores_1ubq_poseop();
	}

	/// @brief test what happens when we do not use FillAUTO (meaning - is there no AUTO?)
	void test_noAUTO() {

		core::pack::task::TaskFactory tf;
		//tf.push_back(core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::FillAUTOTaskOperation() ) );

		core::pack::task::PackerTaskOP task( tf.create_task_and_apply_taskoperations( *poseop ) );

		for ( core::Size i = 1; i <= poseop->size(); ++i ) {
			TS_ASSERT(!task->residue_task(i).has_behavior("AUTO"));
		}
	}

	void test_FillAUTO() {

		core::pack::task::TaskFactory tf;
		tf.push_back(core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::FillAUTOTaskOperation() ) );

		core::pack::task::PackerTaskOP task( tf.create_task_and_apply_taskoperations( *poseop ) );

		for ( core::Size i = 1; i <= poseop->size(); ++i ) {
			TS_ASSERT(task->residue_task(i).has_behavior("AUTO"));
		}
	}




};
