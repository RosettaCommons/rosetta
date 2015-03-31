// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocals/moves/SimulatedTempering.cc
/// @brief Light-weight class for simulated tempering.
/// @author Fang-Chieh Chou

// Unit Headers
#include <protocols/moves/SimulatedTempering.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Numeric Headers
#include <numeric/random/random.hh>

// Utility Headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cmath>

static thread_local basic::Tracer TR( "protocols.moves.SimulatedTempering" );

using namespace core;
using namespace core::pose;
using namespace core::scoring;

namespace protocols {
namespace moves {


SimulatedTempering::SimulatedTempering(
	Pose & pose,
	ScoreFunctionOP const & scorefxn,
	utility::vector1<Real> const & temperatures,
	utility::vector1<Real> const & weights
):
	temperatures_( temperatures ),
	weights_( weights ),
	scorefxn_( scorefxn ),
	rep_scorefxn_( core::scoring::ScoreFunctionOP( new ScoreFunction ) ),
	temp_id_( 1 ),
	cached_score_( ( *scorefxn )( pose ) ),
	rep_cutoff_( 0 )
{
	runtime_assert( temperatures.size() == weights.size() );
	rep_scorefxn_->set_weight( fa_rep, scorefxn->get_weight( fa_rep ) );
}


bool SimulatedTempering::boltzmann( Pose & pose ) {
	Real new_score( 0 );

	if ( rep_cutoff_ > 0 ) {
		// If the repulsion score surpass the cutoff, just skip the full scoring
		// and use repulsion score.
		new_score = ( *rep_scorefxn_ )( pose );
		if ( new_score < rep_cutoff_ ) new_score = ( *scorefxn_ )( pose );
	} else {
		new_score = ( *scorefxn_ )( pose );
	}

	// Hack: negative temperature is infinite T.
	if ( temperature() < 0 ) {
		cached_score_ = new_score;
		return true;
	}

	// Metropolis criterion
	if (
		new_score <= cached_score_ ||
		numeric::random::rg().uniform() < exp( ( cached_score_ - new_score ) / temperature() )
	) {
		cached_score_ = new_score;
		return true;
	}
	return false;
}


bool SimulatedTempering::t_jump() {
	// Early exit if only 1 temperature exists
	if ( temperatures_.size() == 1 ) return false;

	// Obtain the new temperature id. Only move one step upward or downward.
	Size new_temp_id;
	if ( temperatures_.size() == 2 ) {
		// Special case: only 2 temperatures.
		new_temp_id = ( temp_id_ == 1 ) ? 2 : 1;
	} else {
		if ( numeric::random::rg().uniform() < 0.5 ) {
			if ( temp_id_ == temperatures_.size() ) return false;  // Out of bound
			new_temp_id = temp_id_ + 1;
		} else {
			if ( temp_id_ == 1 ) return false; // out of bound
			new_temp_id = temp_id_ - 1;
		}
	}

	// Now check if accept the new id
	Real inv_t_curr = 1.0 / temperature();
	Real inv_t_new = 1.0 / temperatures_[new_temp_id];
	// Convert negative T to inf.
	if ( inv_t_curr < 0 ) inv_t_curr = 0;
	if ( inv_t_new < 0 ) inv_t_new = 0;

	Real const log_prob(
		cached_score_ * inv_t_curr - cached_score_ * inv_t_new +
		weights_[new_temp_id] - weights_[temp_id_]
	);

	if ( log_prob >= 0 || numeric::random::rg().uniform() < exp( log_prob ) ) {
		temp_id_ = new_temp_id;
		return true;
	}
	return false;
}


void SimulatedTempering::score_function( ScoreFunctionOP const & scorefxn ) {
	scorefxn_ = scorefxn;
}


ScoreFunctionOP SimulatedTempering::score_function() const {
	return scorefxn_;
}

} // namespace moves
} // namespace protocols

