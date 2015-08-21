// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/refiners/RotamerTrialsRefiner.hh>
#include <protocols/loop_modeling/refiners/RotamerTrialsRefinerCreator.hh>
#include <protocols/loop_modeling/refiners/packing_helper.hh>
#include <protocols/loop_modeling/utilities/rosetta_scripts.hh>

// Core headers
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// RosettaScripts headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace loop_modeling {
namespace refiners {

using namespace std;
using core::pose::Pose;
using core::pack::task::TaskFactoryOP;
using core::scoring::ScoreFunctionOP;

protocols::moves::MoverOP RotamerTrialsRefinerCreator::create_mover() const {
	return protocols::moves::MoverOP( new RotamerTrialsRefiner );
}

std::string RotamerTrialsRefinerCreator::keyname() const {
	return "RotamerTrialsRefiner";
}

void RotamerTrialsRefiner::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose) {

	LoopMover::parse_my_tag(tag, data, filters, movers, pose);
	utilities::set_scorefxn_from_tag(*this, tag, data);
	utilities::set_task_factory_from_tag(*this, tag, data);
}

bool RotamerTrialsRefiner::do_apply(Pose & pose) {
	return packing_helper(pose, this, core::pack::rotamer_trials);
}

ScoreFunctionOP RotamerTrialsRefiner::get_score_function() {
	return get_tool<ScoreFunctionOP>(ToolboxKeys::SCOREFXN);
}

void RotamerTrialsRefiner::set_score_function(ScoreFunctionOP score_function) {
	set_tool(ToolboxKeys::SCOREFXN, score_function);
}

TaskFactoryOP RotamerTrialsRefiner::get_task_factory() {
	return get_tool<TaskFactoryOP>(ToolboxKeys::TASK_FACTORY);
}

TaskFactoryOP RotamerTrialsRefiner::get_task_factory(TaskFactoryOP fallback) {
	return get_tool<TaskFactoryOP>(ToolboxKeys::TASK_FACTORY, fallback);
}

void RotamerTrialsRefiner::set_task_factory(TaskFactoryOP task_factory) {
	set_tool(ToolboxKeys::TASK_FACTORY, task_factory);
}

}
}
}
