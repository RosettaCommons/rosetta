// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_utilities_PeriodicMover_HH
#define INCLUDED_protocols_loop_modeling_utilities_PeriodicMover_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/utilities/PeriodicMover.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

/// @brief Apply another mover at regular intervals.
class PeriodicMover : public LoopMover {

public:

	/// @brief Constructor with mover and period arguments.
	PeriodicMover(LoopMoverOP mover, Size period);

	/// @copydoc LoopMover::get_name
	string get_name() const { return "PeriodicMover"; }

protected:

	/// @brief Apply the wrapped mover once every `period` iterations.
	bool do_apply(Pose & pose);

private:
	LoopMoverOP mover_;
	core::Size period_;
	core::Size iteration_;

};

}
}
}

#endif
