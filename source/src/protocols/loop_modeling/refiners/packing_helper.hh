// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_modeling_refiners_utilities_HH
#define INCLUDED_protocols_loop_modeling_refiners_utilities_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>

// Core headers
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loops/Loops.hh>
#include <protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.hh>

namespace protocols {
namespace loop_modeling {
namespace refiners {

template <typename PackingRefiner, typename PackingAlgorithm>
bool packing_helper(
	core::pose::Pose & pose,
	PackingRefiner * refiner,
	PackingAlgorithm packer) {

	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using protocols::toolbox::task_operations::RestrictToLoopsAndNeighbors;
	using protocols::toolbox::task_operations::RestrictToLoopsAndNeighborsOP;

	// See if this refiner was given a custom task factory.

	TaskFactoryOP task_factory = refiner->get_task_factory(TaskFactoryOP());

	// If not, create a default one that just packs within 10A of the loops.

	if ( ! task_factory ) {
		RestrictToLoopsAndNeighborsOP restrict_to_loops( new RestrictToLoopsAndNeighbors );
		restrict_to_loops->set_loops(refiner->get_loops());
		restrict_to_loops->set_include_neighbors(true);
		restrict_to_loops->set_cutoff_distance(10);

		ExtraRotamersGenericOP extra_rotamers( new ExtraRotamersGeneric );
		extra_rotamers->ex1(true);
		extra_rotamers->ex2(true);
		extra_rotamers->extrachi_cutoff(0);

		task_factory = TaskFactoryOP( new TaskFactory );
		task_factory->push_back(TaskOperationCOP( new RestrictToRepacking ));
		task_factory->push_back(restrict_to_loops);
		task_factory->push_back(extra_rotamers);
	}

	// Repack the pose.

	PackerTaskOP packer_task = task_factory->create_task_and_apply_taskoperations(pose);
	packer_task->set_bump_check(true);
	packer(pose, *refiner->get_score_function(), packer_task);
	pose.update_residue_neighbors();

	return true;
}

}
}
}

#endif
