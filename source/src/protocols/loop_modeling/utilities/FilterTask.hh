// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_utilities_FilterTask_HH
#define INCLUDED_protocols_loop_modeling_utilities_FilterTask_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMoverTask.hh>
#include <protocols/loop_modeling/utilities/FilterTask.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocols headers
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

class FilterTask : public LoopMoverTask {

public:
	FilterTask(protocols::filters::FilterOP filter);
	string get_name() const { return "FilterTask"; }

public:
	bool apply(Pose & pose, Loop const & loop, ScoreFunctionCOP score_function);

private:
	protocols::filters::FilterOP filter_;

};

}
}
}


#endif

