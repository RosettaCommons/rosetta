#pragma once

// Unit headers
#include <protocols/loop_modeling/LoopMoverTask.hh>
#include <protocols/loop_modeling/filters/BumpCheck.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocols headers
#include <protocols/loops/Loop.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace filters {

using core::pose::Pose;
using core::scoring::ScoreFunctionCOP;
using protocols::loops::Loop;
using protocols::loop_modeling::LoopMoverTask;

class RamaCheck : public LoopMoverTask {

public:
	bool apply(Pose & pose, Loop const & loop, ScoreFunctionCOP score_function);
	string get_name() const { return "RamaCheck"; }

};

}
}
}


