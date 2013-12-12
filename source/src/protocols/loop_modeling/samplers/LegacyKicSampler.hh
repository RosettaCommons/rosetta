// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_samplers_LegacyKicSampler_HH
#define INCLUDED_protocols_loop_modeling_samplers_LegacyKicSampler_HH

// Unit headers
#include <protocols/loop_modeling/LoopMoverTask.hh>
#include <protocols/loop_modeling/samplers/LegacyKicSampler.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace samplers {

using core::pose::Pose;
using core::scoring::ScoreFunctionCOP;
using protocols::loops::Loop;
using protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP;

class LegacyKicSampler : public LoopMoverTask {

public:
	LegacyKicSampler();
	string get_name() const { return "LegacyKicSampler"; }

public:
	bool apply(Pose & pose, Loop const & loop, ScoreFunctionCOP score_function);

private:
	KinematicMoverOP mover_;

};

}
}
}


#endif

