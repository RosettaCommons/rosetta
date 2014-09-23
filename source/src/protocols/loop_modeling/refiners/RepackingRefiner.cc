// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/refiners/RepackingRefiner.hh>
#include <protocols/loop_modeling/refiners/RepackingRefinerCreator.hh>

// Core headers
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_main.hh>

// Utility headers
#include <utility/vector1.hh>
#include <iostream>

namespace protocols {
namespace loop_modeling {
namespace refiners {

using namespace std;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using protocols::loops::Loop;

protocols::moves::MoverOP RepackingRefinerCreator::create_mover() const {
	return new RepackingRefiner;
}

std::string RepackingRefinerCreator::keyname() const {
	return "RepackingRefiner";
}

RepackingRefiner::RepackingRefiner() {
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using core::pack::task::operation::InitializeFromCommandline;
	using core::pack::task::operation::IncludeCurrent;

	task_factory_ = new TaskFactory;
	task_factory_->push_back(new InitializeFromCommandline);
	task_factory_->push_back(new IncludeCurrent);

	// Should read resfile if present.
}

bool RepackingRefiner::do_apply(Pose & pose) {
	using core::pack::pack_rotamers;
	using core::pack::task::PackerTaskOP;

	PackerTaskOP packer_task =
		task_factory_->create_task_and_apply_taskoperations(pose);

	utility::vector1<bool> loop_residues = 
		protocols::loops::select_loop_residues(pose, get_loops(), true, 10.0);

	core::pose::symmetry::make_residue_mask_symmetric(pose, loop_residues);

	packer_task->set_bump_check(true);
	packer_task->restrict_to_repacking();
	packer_task->restrict_to_residues(loop_residues);

	pack_rotamers(pose, *get_score_function(), packer_task);
	pose.update_residue_neighbors();

	return true;
}

}
}
}


