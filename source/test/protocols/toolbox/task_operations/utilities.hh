// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// C++ headers
#include <string>
#include <iostream>


inline
void
test_task_operation(
		core::pack::task::operation::TaskOperationCOP task_op,
		core::pose::Pose pose,
		std::string packable_mask,
		std::string designable_mask) {

	using core::Size;

	core::pack::task::TaskFactory task_factory; task_factory.push_back(task_op);
	core::pack::task::PackerTaskOP task =
		task_factory.create_task_and_apply_taskoperations(pose);

	for (Size i = 1; i <= pose.total_residue(); i++) {
		bool expecting_to_pack = (packable_mask[i-1] == '1');
		TS_ASSERT_EQUALS(expecting_to_pack, task->being_packed(i));
	}

	for (Size i = 1; i <= pose.total_residue(); i++) {
		bool expecting_to_design = (designable_mask[i-1] == '1');
		TS_ASSERT_EQUALS(expecting_to_design, task->being_designed(i));
	}
}


