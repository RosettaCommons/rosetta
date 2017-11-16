// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/canonical_sampling/TemperatureControllerMover.cc
/// @brief TemperatureController methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/TemperatureController.hh>


// protocols headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>

#include <protocols/rosetta_scripts/util.hh>

// core headers
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/file/file_sys_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cmath>

static basic::Tracer tr( "protocols.canonical_sampling.TemperatureController" );

namespace protocols {
namespace canonical_sampling {
using namespace core;


std::string
interpolation_type_enum_to_string( InterpolationType interp_enum ) {

	return interp_enum == linear ? "linear" : "exponential";
}

InterpolationType
interpolation_type_string_to_enum( std::string const & interp_string ) {

	if ( interp_string == "linear" ) {
		return linear;
	} else if ( interp_string != "exponential" ) {
		utility_exit_with_message("invalid temp_interpolation value, expecting linear or exponential");
	}
	return exponential;
}

TemperatureController::TemperatureController() :
	protocols::canonical_sampling::ThermodynamicObserver()
{}

TemperatureController::TemperatureController( TemperatureController const & other ) :
	protocols::canonical_sampling::ThermodynamicObserver(other)
{}

void
TemperatureController::observe_after_metropolis( protocols::canonical_sampling::MetropolisHastingsMover const& ) {
	Real const score( monte_carlo_->last_accepted_score() );
	temperature_move( score );
}

protocols::moves::MonteCarloCOP
TemperatureController::monte_carlo() const {
	return monte_carlo_;
}

void
TemperatureController::set_monte_carlo(
	protocols::moves::MonteCarloOP monte_carlo
)
{
	monte_carlo_ = monte_carlo;
}

protocols::moves::MonteCarloOP
TemperatureController::monte_carlo() {
	return monte_carlo_;
}

core::Real
TemperatureController::temperature_move( core::pose::Pose& pose ) {
	return temperature_move( pose.energies().total_energy() );
}

std::string
TemperatureController::get_name() const
{
	return "TemperatureController";
}


void
TemperatureController::initialize_simulation(
	core::pose::Pose &,
	protocols::canonical_sampling::MetropolisHastingsMover const &,
	core::Size level,
	core::Real,
	core::Size
) {
	if ( level != 1 ) {
		utility_exit_with_message( "TemperatureController - Baseclass doesn't know about different temp-levels" );
	};
}


} //moves
} //protocols

