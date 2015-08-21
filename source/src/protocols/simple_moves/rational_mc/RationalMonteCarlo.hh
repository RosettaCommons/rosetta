// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/rational_mc/RationalMonteCarlo.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_rational_mc_RationalMonteCarlo_HH
#define INCLUDED_protocols_simple_moves_rational_mc_RationalMonteCarlo_HH

// Unit header
#include <protocols/simple_moves/rational_mc/RationalMonteCarlo.fwd.hh>

// C/C++ headers
#include <string>

// External headers
#include <boost/function.hpp>
#include <boost/unordered/unordered_map.hpp>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Package headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace simple_moves {
namespace rational_mc {

/// @brief Trigger API definition
typedef boost::function<void(const core::pose::Pose&)> RationalMonteCarloTrigger;
typedef boost::unordered_map<core::Size, RationalMonteCarloTrigger> Triggers;

/// @class Trial-based Monte Carlo minization primitive. Do not modify this
/// class; almost anything you could possibly want to add is a bad idea.
class RationalMonteCarlo : public protocols::moves::Mover {
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;

public:
	RationalMonteCarlo(moves::MoverOP mover,
		ScoreFunctionOP score,
		core::Size num_trials,
		core::Real temperature,
		bool recover_low);

	/// @brief Applies the underlying mover to <pose> the specified number of times
	void apply(Pose& pose);

	/// @brief Returns this mover's name
	std::string get_name() const;

	// -- Accessors -- //
	bool recover_low() const;
	moves::MoverOP mover() const;
	core::Size num_trials() const;
	core::Real temperature() const;
	const core::scoring::ScoreFunction& score_function() const;

	/// @brief Updates the last accepted and lowest scoring pose members of the
	/// MonteCarlo member variable to the score of the specified pose
	void reset(const core::pose::Pose& pose);
	const core::pose::Pose& lowest_score_pose() const;
	const core::pose::Pose& last_accepted_pose() const;

	// -- Mutators -- //
	void set_mover(moves::MoverOP mover);
	void set_recover_low(bool recover_low);
	void set_num_trials(core::Size num_trials);
	void set_temperature(core::Real temperature);

	void enable_autotemp(core::Real quench);
	void disable_autotemp();

	/// @brief Updates the score function. Before calling this method, make sure
	/// that apply() has been called at least once on a non-empty pose.
	void set_score_function(core::scoring::ScoreFunctionOP score);

	// -- Triggers -- //

	/// @brief Registers the specified trigger with this instance. Returns a unique
	/// identifier for referring to this trigger in subsequent operations.
	core::Size add_trigger(const RationalMonteCarloTrigger& trigger);

	/// @brief Unregisters the trigger with the given unique identifier
	void remove_trigger(core::Size trigger_id);

protected:
	/// @brief Executes all triggers attached to this instance. The order of trigger
	/// execution is undefined. Do not assume, depend, or in any way rely upon a
	/// partiular ordering.
	void fire_all_triggers(const Pose& pose);

private:
	/// @brief Underlying mover
	moves::MoverOP mover_;

	/// @brief Determines whether applications of the base mover should be accepted
	/// by applying the Metropolis criterion to the pose
	moves::MonteCarloOP mc_;

	/// @brief Number of times to execute the underlying mover in calls to apply()
	core::Size num_trials_;

	/// @brief Determines whether the low scoring pose should be recovered at the
	/// conclusion of the apply() method
	bool recover_low_;

	/// @brief Collection of function callbacks
	Triggers triggers_;

	/// @brief Next trigger id to be assigned
	core::Size next_trigger_id_;
};

}  // namespace rational_mc
}  // namespace simple_moves
}  // namespace protocols

#endif  // protocols_simple_moves_rational_mc_RationalMonteCarlo_HH_
