// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_utilities_RepeatedMover_HH
#define INCLUDED_protocols_loop_modeling_utilities_RepeatedMover_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/utilities/RepeatedMover.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

/// @brief Apply another mover multiple times.
class RepeatedMover : public LoopMover {

public:

	/// @brief Constructor with mover and repeat arguments.
	RepeatedMover(LoopMoverOP mover, Size iterations);

	/// @copydoc LoopMover::get_name
	string get_name() const { return "RepeatedMover"; }

protected:

	/// @brief Apply the wrapped move `iterations` times.
	bool do_apply(Pose & pose);

private:
	LoopMoverOP mover_;
	Size iterations_;

};

}
}
}

#endif
