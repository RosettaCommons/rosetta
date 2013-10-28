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
#include <protocols/canonical_sampling/TrialCounterObserverCreator.hh>

#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
// Package Headers

// Project Headers
#include <core/pose/Pose.hh>


// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

// C++ Headers

static basic::Tracer tr( "protocols.canonical_sampling.TrialCounter" );

namespace protocols {
namespace canonical_sampling {


std::string
TrialCounterObserverCreator::keyname() const {
	return TrialCounterObserverCreator::mover_name();
}

protocols::moves::MoverOP
TrialCounterObserverCreator::create_mover() const {
	return new TrialCounterObserver;
}

std::string
TrialCounterObserverCreator::mover_name() {
	return "TrialCounterObserver";
}


///@brief
TrialCounterObserver::TrialCounterObserver(
) : ThermodynamicObserver()
{
	Mover::type( "TrialCounterObserver" );
}

TrialCounterObserver::~TrialCounterObserver() {}

std::string TrialCounterObserver::get_name() const {
	return "TrialCounterObserver";
}

protocols::moves::MoverOP
TrialCounterObserver::clone() const {
	return new TrialCounterObserver( *this );
}



void
TrialCounterObserver::parse_my_tag(
       utility::tag::TagCOP const tag,
       basic::datacache::DataMap &,
       protocols::filters::Filters_map const &,
       protocols::moves::Movers_map const &,
       core::pose::Pose const &
) { //no options ...
	file_ = tag->getOption< std::string >("file","trial.stats");
}


void
TrialCounterObserver::initialize_simulation(
	core::pose::Pose & /*pose*/,
	MetropolisHastingsMover const& mhm, /*metropolis_hastings_mover*/
	core::Size //default=0; non-zero if trajectory is restarted
)
{
	counters_.set_temperature_observer( mhm.tempering() );
	counters_.reset();
}

	/// @brief callback executed after the Metropolis criterion is evaluated
void
TrialCounterObserver::observe_after_metropolis(
		MetropolisHastingsMover const & metropolis_hastings_mover
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
	MetropolisHastingsMover const &mhm /*metropolis_hastings_mover*/
)
{
	counters_.show( tr.Info );
	counters_.write_to_file( file_, mhm.output_name() );
}

} //moves
} //protocols

