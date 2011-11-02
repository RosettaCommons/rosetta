// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/moves/ThermodynamicObserver.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_moves_ThermodynamicObserver_hh
#define INCLUDED_protocols_moves_ThermodynamicObserver_hh

// Unit Headers
#include <protocols/moves/ThermodynamicObserver.fwd.hh>

// Project Headers
#include <protocols/moves/MetropolisHastingsMover.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Range.hh>
// AUTO-REMOVED #include <core/id/TorsionID_Range.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

///@details
class ThermodynamicObserver : public protocols::moves::Mover {

public:

	///@brief
	ThermodynamicObserver();

	virtual
	~ThermodynamicObserver();

	/// @brief callback executed before any Monte Carlo trials
	virtual
	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	/// @brief callback executed after the Metropolis criterion is evaluated
	virtual
	void
	observe_after_metropolis(
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
	) = 0;

	/// @brief callback executed after all Monte Carlo trials
	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
	);

private:

}; //end ThermodynamicObserver

} //namespace moves
} //namespace protocols

#endif // INCLUDED_protocols_moves_ThermodynamicObserver_HH
