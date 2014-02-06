// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopMoverTask_HH
#define INCLUDED_protocols_loop_modeling_LoopMoverTask_HH

// Unit headers
#include <protocols/loop_modeling/LoopMoverTask.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <boost/utility.hpp>

// C++ headers
#include <string>

namespace protocols {
namespace loop_modeling {

using namespace std;
using utility::pointer::ReferenceCount;
using boost::noncopyable;

class LoopMoverTask : public ReferenceCount, protected noncopyable {

public: 

	/// @brief Return the name of this task.
	virtual string get_name() const = 0;
	
	/// @brief Setup the task.
	virtual void setup(
			core::pose::Pose &,
			protocols::loops::Loop const &,
			core::scoring::ScoreFunctionOP) {}

	/// @brief Make a change to the given pose and/or decide whether or not 
	/// further processing should be done.
	virtual bool apply(
			core::pose::Pose & pose,
			protocols::loops::Loop const & loop,
			core::scoring::ScoreFunctionCOP score_function) = 0;

	/// @brief Allow the task to respond 
	virtual void debrief(bool /*was_accepted*/) {}

};

/// @brief Easy way to make a repeated task.
LoopMoverTaskOP operator * (LoopMoverTaskOP task, core::Size iterations);

/// @brief Easy way to make a periodic task.
LoopMoverTaskOP operator % (LoopMoverTaskOP task, core::Size period);

}
}


#endif

