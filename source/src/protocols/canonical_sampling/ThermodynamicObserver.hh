// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/canonical_sampling/ThermodynamicObserver.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_ThermodynamicObserver_hh
#define INCLUDED_protocols_canonical_sampling_ThermodynamicObserver_hh

// Unit Headers
#include <protocols/canonical_sampling/ThermodynamicObserver.fwd.hh>

// Project Headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace canonical_sampling {

/// @brief Base class for reporting and recording data from a simulation.
class ThermodynamicObserver : public protocols::moves::Mover {

public:

	/// @brief Default constructor.
	ThermodynamicObserver();

	/// @brief Destructor.
	virtual
	~ThermodynamicObserver();

	/// @brief Callback executed after each move is made.
	/// @details Even though the argument is a reference to a non-const pose, 
	/// this method should not make any changes to the pose.  Making changes to 
	/// the pose is the role of the ThermodynamicMover class.  The role of this 
	/// class is to simply observe the poses being generated.
	virtual
	void apply( core::pose::Pose& ) {};

	/// @brief Callback executed before any Monte Carlo trials are attempted.
	virtual
	void
	initialize_simulation(
		core::pose::Pose &,
		MetropolisHastingsMover const &,
		core::Size //non-zero if trajectory is restarted
	) {};

	/// @brief Callback executed after the Metropolis criterion is evaluated.
	virtual
	void
	observe_after_metropolis(
		MetropolisHastingsMover const & metropolis_hastings_mover
	) = 0;

	/// @brief Callback executed after all Monte Carlo trials are completed.
	virtual
	void
	finalize_simulation(
		core::pose::Pose &,
		MetropolisHastingsMover const &
	) {};

	/// @brief Attempt to restart the last simulation that was recorded by this 
	/// observer.
	///
	/// @details For example, consider an observer that records trajectories.  
	/// This method should open the file that was going to be written, read out 
	/// the last pose in that trajectory, and assign it to the given pose 
	/// reference so that the current trajectory can start from the same place.  
	/// Other observers may help setup other parts of the simulation.
	///
	/// This is not a particularly robust system, because it may require several 
	/// unrelated observers working in concert to properly reconstitute the 
	/// simulation.  In fact, the restart feature in MetropolisHastingsMover is 
	/// currently commented out, so this method is never actually invoked.  I 
	/// would advise reimplementing this method to utility_exit_with_message() in 
	/// any subclasses you write, so that you don't waste time writing an unused 
	/// method but so you don't confuse anyone if this feature gets revived in 
	/// the future.
	virtual
	bool
	restart_simulation(
		core::pose::Pose &,
		MetropolisHastingsMover&,
		core::Size&,
		core::Size&,
		core::Real& 
	) { return false; }

	/// @brief Return false if this observer does not require a valid pose.  
	/// TrialCounterObserver is an example of such an observer.
	virtual
	bool
	requires_pose() { return true; };

private:

}; //end ThermodynamicObserver

} //namespace canonical_sampling
} //namespace protocols

#endif // INCLUDED_protocols_canonical_sampling_ThermodynamicObserver_HH
