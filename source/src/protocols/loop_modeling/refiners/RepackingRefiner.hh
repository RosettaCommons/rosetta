// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_refiners_RepackingRefiner_HH
#define INCLUDED_protocols_loop_modeling_refiners_RepackingRefiner_HH

// Unit headers
#include <protocols/loop_modeling/LoopMoverTask.hh>
#include <protocols/loop_modeling/refiners/RepackingRefiner.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocols headers
#include <protocols/loops/Loop.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace refiners {

using core::pose::Pose;
using core::pack::task::TaskFactoryOP;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using protocols::loops::Loop;

class RepackingRefiner : public LoopMoverTask {

public:

	RepackingRefiner();
	string get_name() const { return "RepackingRefiner"; }

	void setup(Pose & pose, Loop const & loop, ScoreFunctionOP score_function);
	bool apply(Pose & pose, Loop const & loop, ScoreFunctionCOP score_function);

private:

	TaskFactoryOP task_factory_;

};

}
}
}


#endif

