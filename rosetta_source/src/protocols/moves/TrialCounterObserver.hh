// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/moves/TrialCounterObserver.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_moves_TrialCounterObserver_hh
#define INCLUDED_protocols_moves_TrialCounterObserver_hh

// Unit Headers
//#include <protocols/moves/TrialCounterObserver.fwd.hh>

// Project Headers
#include <protocols/moves/MetropolisHastingsMover.fwd.hh>
#include <protocols/moves/ThermodynamicObserver.hh>
#include <protocols/moves/MultiTemperatureTrialCounter.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>

namespace protocols {
namespace moves {

///@details
class TrialCounterObserver : public ThermodynamicObserver {

public:

	///@brief
	TrialCounterObserver();

	virtual
	~TrialCounterObserver();

	virtual std::string get_name() const;
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
	);

	/// @brief callback executed after all Monte Carlo trials
	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	/// @brief return false here if a valid pose is not required for "observe"
	/// i.e. a trialcounter
	virtual
	bool
	requires_pose() { return false; }

private:
	MultiTemperatureTrialCounter counters_;
}; //end TrialCounterObserver

} //namespace moves
} //namespace protocols

#endif // INCLUDED_protocols_moves_TrialCounterObserver_HH
