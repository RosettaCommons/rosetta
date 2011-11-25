// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/TrialCounterObserver.cc
/// @brief TrialCounterObserver methods implemented
/// @author


// Unit Headers
#include <protocols/moves/TrialCounterObserver.hh>

#include <protocols/moves/MetropolisHastingsMover.hh>
// Package Headers

// Project Headers
#include <core/pose/Pose.hh>



// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/exit.hh>

// C++ Headers

static basic::Tracer tr( "protocols.moves.TrialCounter" );

namespace protocols {
namespace moves {

///@brief
TrialCounterObserver::TrialCounterObserver(
) : ThermodynamicObserver()
{
	Mover::type( "TrialCounterObserver" );
}

TrialCounterObserver::~TrialCounterObserver() {}

std::string TrialCounterObserver::get_name() const {
	return "TrialCounterObserver";
};

void
TrialCounterObserver::initialize_simulation(
	core::pose::Pose & /*pose*/,
	protocols::moves::MetropolisHastingsMover const& mhm /*metropolis_hastings_mover*/
)
{
	counters_.set_temperature_observer( mhm.tempering() );
	counters_.reset();
}

	/// @brief callback executed after the Metropolis criterion is evaluated
void
TrialCounterObserver::observe_after_metropolis(
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
) {
	std::string const& move_type( metropolis_hastings_mover.last_move().type() );
	counters_.count_trial( move_type );
	if ( metropolis_hastings_mover.last_accepted() ) {
		counters_.count_accepted( move_type );
	}
}

void
TrialCounterObserver::finalize_simulation(
	core::pose::Pose & /*pose*/,
	protocols::moves::MetropolisHastingsMover const &mhm /*metropolis_hastings_mover*/
)
{
	counters_.show( tr.Info );
	counters_.write_to_file( "trial_counts.stats", mhm.output_name() );
}

} //moves
} //protocols

