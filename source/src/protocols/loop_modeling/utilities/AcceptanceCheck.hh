// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_utilities_AcceptanceCheck_HH
#define INCLUDED_protocols_loop_modeling_utilities_AcceptanceCheck_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/utilities/AcceptanceCheck.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Protocols headers
#include <protocols/moves/MonteCarlo.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

/// @brief Apply an additional Monte Carlo acceptance check.
///
/// @details This class is meant to be used when you have a complex loop
/// sampling protocol (e.g. when you are stringing together a large number of
/// loop movers) and you want to add intermediate acceptance checks.  Note that
/// you don't usually have to use this class, because LoopProtocol always makes
/// an acceptance check after all the loop movers have been applied.

class AcceptanceCheck : public LoopMover {

public:
	/// @brief Constructor with optional name argument.
	AcceptanceCheck(
		protocols::moves::MonteCarloOP monte_carlo,
		string name="loop_move");

	/// @copydoc LoopMover::get_name
	string get_name() const { return "AcceptanceCheck"; }

protected:
	/// @brief Apply the Metropolis criterion to the given pose and return true
	/// if the pose was accepted.
	bool do_apply(Pose & pose);

private:
	protocols::moves::MonteCarloOP monte_carlo_;
	string name_;

};

}
}
}


#endif

