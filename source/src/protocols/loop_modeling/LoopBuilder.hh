// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopBuilder_HH
#define INCLUDED_protocols_loop_modeling_LoopBuilder_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/LoopBuilder.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocol headers
#include <protocols/kinematic_closure/KicMover.fwd.hh>
#include <protocols/loop_modeling/refiners/LocalMinimizationRefiner.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>

namespace protocols {
namespace loop_modeling {

/// @brief Build loops from scratch.
///
/// @details Building a loop from scratch is useful in two scenarios.  The 
/// first is when there's missing density that needs to be modeled, and the 
/// second is when the whole loop modeling algorithm needs to be benchmarked.  
/// To build a loop, every torsion in the loop is picked randomly from the 
/// Ramachandran distribution and every other DOF is set to its ideal value.  
/// Kinematic closure is then used to connect the backbone.  If the resulting 
/// structure passes a rather lenient bump check, building is complete.  If 
/// not, the structure is thrown out and the whole process is tried again.
///
/// This process can be very slow for long loops, because there's nothing 
/// guiding the algorithm towards the right solution.  Torsions are just being 
/// randomly picked, and very often they won't fit in the relatively narrow 
/// space that's available.  The problem is worse for interior loops than it is 
/// for surface loops, of course.   This algorithm seems to work well enough on 
/// 12 residue loops, but beyond that it may be necessary to develop a smarter 
/// algorithm preferentially builds into free space.

class LoopBuilder : public LoopMover {

public:

	/// @brief Default constructor.
	LoopBuilder();

	/// @copydoc LoopMover::get_name
	string get_name() const { return "LoopBuilder"; }

protected:

	/// @brief Attempt to find a reasonable loop conformation without using any 
	/// information from the original coordinates.
	bool do_apply(Pose & pose, Loop const & loop);

private:

	protocols::kinematic_closure::KicMoverOP kic_mover_;
	refiners::LocalMinimizationRefinerOP minimizer_;

};

}
}

#endif

