// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_kinematic_closure_pivot_pickers_FixedOffsetsPivots_HH
#define INCLUDED_protocols_kinematic_closure_pivot_pickers_FixedOffsetsPivots_HH

// Unit headers
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.hh>
#include <protocols/kinematic_closure/pivot_pickers/FixedOffsetsPivots.fwd.hh>

namespace protocols {
namespace kinematic_closure {
namespace pivot_pickers {

/// @brief Randomly pick pivots which are always offsets given in a list but
/// which can appear anywhere in the loop.
class FixedOffsetsPivots : public PivotPicker {

public:
	/// @brief Default constructor.
	FixedOffsetsPivots(utility::vector1<Size> offsets) : offsets_(offsets) {}

public:
	/// @copydoc PivotPicker::get_name
	std::string get_name() const { return "FixedOffsetsPivots"; }

	/// @copydoc PivotPicker::pick
	Loop pick(Pose const & pose, Loop const & loop);

private:
	utility::vector1<Size> offsets_;

};

}
}
}

#endif

