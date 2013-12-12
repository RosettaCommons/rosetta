// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_kinematic_closure_pivot_pickers_FixedPivots_HH
#define INCLUDED_protocols_kinematic_closure_pivot_pickers_FixedPivots_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.hh>
#include <protocols/kinematic_closure/pivot_pickers/FixedPivots.fwd.hh>

namespace protocols {
namespace kinematic_closure {
namespace pivot_pickers {

/// @brief Use a fixed set of pivots specified in advance.
/// @details The pivots are specified in the constructor.  Make sure that the 
/// pivots will always be contained in the loop being sampled, otherwise 
/// strange behavior may occur.

class FixedPivots : public PivotPicker {

public:
	/// @brief Constructor which specifies the fixed pivots.
	FixedPivots(Size start, Size stop, Size cut);

public:
	/// @copydoc PivotPicker::pick
	Loop pick(Pose const & pose, Loop const & loop);

private:
	Size start_, stop_, cut_;

};

}
}
}

#endif

