// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/canonical_sampling/TrialCounterObserver.cc
/// @brief protocols::canonical_sampling::TrialCounterObserver methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/TrialCounterObserver.hh>

#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
// Package Headers

// Project Headers
#include <core/pose/Pose.hh>



// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/exit.hh>

// C++ Headers

static basic::Tracer tr( "protocols.canonical_sampling.TrialCounter" );

namespace protocols {
namespace canonical_sampling {

///@brief
protocols::canonical_sampling::TrialCounterObserver::TrialCounterObserver(
) : protocols::canonical_sampling::ThermodynamicObserver()
{
	Mover::type( "TrialCounterObserver" );
}

protocols::canonical_sampling::TrialCounterObserver::~TrialCounterObserver() {}

std::string protocols::canonical_sampling::TrialCounterObserver::get_name() const {
	return "TrialCounterObserver";
};

void
protocols::canonical_sampling::TrialCounterObserver::initialize_simulation(
	core::pose::Pose & /*pose*/,
	protocols::canonical_sampling::MetropolisHastingsMover const& mhm /*metropolis_hastings_mover*/
)
{
	counters_.set_temperature_observer( mhm.tempering() );
	counters_.reset();
}

	/// @brief callback executed after the Metropolis criterion is evaluated
void
protocols::canonical_sampling::TrialCounterObserver::observe_after_metropolis(
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
) {
	std::string const& move_type( metropolis_hastings_mover.last_move().type() );
	counters_.count_trial( move_type );
	if ( metropolis_hastings_mover.last_accepted() ) {
		counters_.count_accepted( move_type );
	}
}

void
protocols::canonical_sampling::TrialCounterObserver::finalize_simulation(
	core::pose::Pose & /*pose*/,
	protocols::canonical_sampling::MetropolisHastingsMover const &mhm /*metropolis_hastings_mover*/
)
{
	counters_.show( tr.Info );
	counters_.write_to_file( "trial_counts.stats", mhm.output_name() );
}

} //moves
} //protocols

