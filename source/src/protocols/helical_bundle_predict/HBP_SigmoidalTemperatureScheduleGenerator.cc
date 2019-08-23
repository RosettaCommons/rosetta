// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBP_SigmoidalTemperatureScheduleGenerator.cc
/// @brief A class to generate a temperature ramping scheme for a simulated annealing trajectory.  This version
/// ramps sigmoidally.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#include <protocols/helical_bundle_predict/HBP_SigmoidalTemperatureScheduleGenerator.hh>

// Basic includes
#include <basic/Tracer.hh>

// Utility includes
#include <utility/pointer/memory.hh>

// STL includes
#include <cmath>

static basic::Tracer TR( "protocols.helical_bundle_predict.HBP_SigmoidalTemperatureScheduleGenerator" );


namespace protocols {
namespace helical_bundle_predict {

/// @brief Default constructor:
HBP_SigmoidalTemperatureScheduleGenerator::HBP_SigmoidalTemperatureScheduleGenerator():
	HBP_TemperatureScheduleGenerator(),
	exponent_(2.0)
{}

/// @brief Options constructor.
HBP_SigmoidalTemperatureScheduleGenerator::HBP_SigmoidalTemperatureScheduleGenerator( core::Real const & max_temp, core::Real const & min_temp ) :
	HBP_TemperatureScheduleGenerator(max_temp, min_temp),
	exponent_(2.0)
{}

/// @brief Copy constructor.
HBP_SigmoidalTemperatureScheduleGenerator::HBP_SigmoidalTemperatureScheduleGenerator( HBP_SigmoidalTemperatureScheduleGenerator const &src ) :
	HBP_TemperatureScheduleGenerator(src),
	exponent_(src.exponent_)
{}

/// @brief Destructor.
HBP_SigmoidalTemperatureScheduleGenerator::~HBP_SigmoidalTemperatureScheduleGenerator(){}

/// @brief Copy this object and return an owning pointe to the copy.
HBP_TemperatureScheduleGeneratorOP
HBP_SigmoidalTemperatureScheduleGenerator::clone() const {
	return HBP_TemperatureScheduleGeneratorOP( utility::pointer::make_shared< HBP_SigmoidalTemperatureScheduleGenerator >( *this ) );
}

///////////////////////////////////////// PUBLIC MEMBER FUNCTIONS /////////////////////////////////////////

/// @brief Calculate the current temperature for this point in the trajectory.
core::Real
HBP_SigmoidalTemperatureScheduleGenerator::calculate_current_temperature(
	core::Size const current_step,
	core::Size const max_steps
) const {
	core::Real const max_temp( get_max_temp_for_current_round() );
	core::Real const min_temp( get_min_temp_for_current_round() );

	debug_assert( max_temp > 0 && min_temp > 0 );
	core::Real const ln_max_temp( std::log( max_temp ) );
	core::Real const ln_min_temp( std::log( min_temp ) );

	core::Real const fract( static_cast<core::Real>( current_step - 1 ) / static_cast<core::Real>( max_steps - 1 ) );

	core::Real const sigmoid(
		fract > 1e-12 ?
		1.0/(1.0 + ( std::pow((1.0-fract)/fract, exponent_) ) ) :
		0.0
	);

	core::Real const tempval( std::exp(ln_min_temp*sigmoid + ln_max_temp*(1.0-sigmoid)) );

	TR << "Fract=" << fract << "\tSigmoid=" << sigmoid << "\tTemperature=" << tempval << std::endl;
	return tempval;
}

/// @brief Set the exponent.
/// @details This is a value that determines how S-shaped the sigmoid is.  A value of 1 gives a
/// straight line; higher values have more lag at start and end.
/// @note Must be greater than zero.
void
HBP_SigmoidalTemperatureScheduleGenerator::set_exponent(
	core::Real const & setting
) {
	runtime_assert_string_msg( setting > 0.0, "Error in HBP_SigmoidalTemperatureScheduleGenerator::set_exponent(): The exponent must be greater than zero." );
	exponent_ = setting;
}

} //protocols
} //helical_bundle_predict






