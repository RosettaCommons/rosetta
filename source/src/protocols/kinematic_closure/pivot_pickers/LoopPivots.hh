// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_kinematic_closure_pivot_pickers_LoopPivots_HH
#define INCLUDED_protocols_kinematic_closure_pivot_pickers_LoopPivots_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.hh>
#include <protocols/kinematic_closure/pivot_pickers/LoopPivots.fwd.hh>

namespace protocols {
namespace kinematic_closure {
namespace pivot_pickers {

/// @brief Use the start, stop, and cut points of the given loop as pivots.
/// @details This algorithm is meant to be simple and intuitive.  The pick() 
/// method will be passed a loop, and the pivots will be taken from the 
/// parameters of that loop.  If the cut point is outside the loop (e.g. if it 
/// hasn't been set), the midpoint of the loop will be used instead.
class LoopPivots : public PivotPicker {

public:
	/// @copydoc PivotPicker::pick
	Loop pick(Pose const & pose, Loop const & loop);

};

}
}
}

#endif

