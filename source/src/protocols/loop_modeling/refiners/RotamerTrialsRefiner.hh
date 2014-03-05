// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_refiners_RotamerTrialsRefiner_HH
#define INCLUDED_protocols_loop_modeling_refiners_RotamerTrialsRefiner_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/refiners/RotamerTrialsRefiner.fwd.hh>

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
using core::scoring::ScoreFunctionCOP;
using protocols::loops::Loop;

/// @brief Refine sampled loops using rotamer trials.
class RotamerTrialsRefiner : public LoopMover {

public:
	
	/// @brief Default constructor.
	RotamerTrialsRefiner();

	/// @copydoc LoopMover::get_name
	string get_name() const { return "RotamerTrialsRefiner"; }

protected:

	/// @brief Use rotamer trials to refine the pose within 10A of the loops 
	/// being sampled.
	bool do_apply(Pose & pose);

private:
	TaskFactoryOP task_factory_;

};

}
}
}


#endif

