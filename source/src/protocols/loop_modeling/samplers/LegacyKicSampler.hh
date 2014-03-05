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
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/samplers/LegacyKicSampler.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace samplers {

/// @brief Apply the legacy KinematicMover.
///
/// @details The old KinematicMover works just as well as the new KicMover, but 
/// it is not written as well and is therefore much harder to understand and 
/// maintain.  Largely for this reason, development will shift towards the new 
/// implementation and the old one will grow more and more out of date.  This 
/// sampler was written to compare the two algorithms and is only meant to be a 
/// debugging tool.  However, it may be useful in the short term while features 
/// like next-gen KIC are still being ported to the new implementation.

class LegacyKicSampler : public LoopMover {

public:
	/// @brief Default constructor.
	LegacyKicSampler();

	/// @copydoc LoopMover::get_name
	string get_name() const { return "LegacyKicSampler"; }

protected:
	/// @brief Apply the legacy KinematicMover algorithm.
	bool do_apply(Pose & pose, Loop const & loop);

private:
	protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP mover_;

};

}
}
}


#endif

