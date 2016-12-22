// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocals/moves/SimulatedTempering.hh
/// @brief Light-weight class for simulated tempering.
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_moves_SimulatedTempering_hh
#define INCLUDED_protocols_moves_SimulatedTempering_hh

// type headers
#include <core/types.hh>

// unit headers
#include <protocols/moves/SimulatedTempering.fwd.hh>

// package headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace protocols {
namespace moves {

class SimulatedTempering : public utility::pointer::ReferenceCount {
public:
	SimulatedTempering(
		core::pose::Pose & pose,
		core::scoring::ScoreFunctionCOP scorefxn,
		utility::vector1<core::Real> const & temperatures,
		utility::vector1<core::Real> const & weights
	);

	/// @brief Applies the Metropolis Criterion on pose.
	bool boltzmann( core::pose::Pose & pose );

	/// @brief Attempt temperature jumping.
	bool t_jump();

	/// @brief Get the id of current temperature
	core::Size temp_id() const { return temp_id_; }

	/// @brief Get the current temperature
	core::Real temperature() const { return temperatures_[temp_id_]; }

	/// @brief Sets the ScoreFunction to  <scorefxn>
	void score_function( core::scoring::ScoreFunctionCOP scorefxn );

	/// @brief Sets cutoff of repusion filter, use 0 to turn it off
	void set_rep_cutoff( core::Real const setting ) { rep_cutoff_ = setting; }

	/// @brief Sets the next boltzmann() call to automatically reject
	void force_next_move_reject() { force_next_move_reject_ = true; }

	/// @brief Returns the MonteCarlo ScoreFunction
	core::scoring::ScoreFunctionCOP score_function() const;

private:
	utility::vector1<core::Real> const temperatures_, weights_;
	core::scoring::ScoreFunctionCOP scorefxn_, rep_scorefxn_;
	core::Size temp_id_;
	core::Real cached_score_, rep_cutoff_;
	bool force_next_move_reject_;

};

} // moves
} // rosetta

#endif
