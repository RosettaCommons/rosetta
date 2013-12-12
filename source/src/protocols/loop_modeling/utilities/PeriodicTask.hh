// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_utilities_PeriodicTask_HH
#define INCLUDED_protocols_loop_modeling_utilities_PeriodicTask_HH

// Unit headers
#include <protocols/loop_modeling/LoopMoverTask.hh>
#include <protocols/loop_modeling/utilities/PeriodicTask.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocols headers
#include <protocols/loops/Loop.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using protocols::loops::Loop;

class PeriodicTask : public LoopMoverTask {

public:

	PeriodicTask(LoopMoverTaskOP, Size period);

	string get_name() const { return "PeriodicTask"; }

public:
	void setup(Pose & pose, Loop const & loop, ScoreFunctionOP score_function);

	bool apply(Pose & pose, Loop const & loop, ScoreFunctionCOP score_function);

	void debrief(bool was_accepted);

private:

	LoopMoverTaskOP task_;
	core::Size iteration_;
	core::Size period_;

};

}
}
}


#endif

