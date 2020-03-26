// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBP_TemperatureScheduleGenerator.cc
/// @brief A base class to generate a temperature ramping scheme for a simulated annealing trajectory.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#include <protocols/helical_bundle_predict/HBP_TemperatureScheduleGenerator.hh>

// Basic includes
#include <basic/Tracer.hh>

// Utility includes
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.helical_bundle_predict.HBP_TemperatureScheduleGenerator" );


namespace protocols {
namespace helical_bundle_predict {

/// @brief Default constructor:
HBP_TemperatureScheduleGenerator::HBP_TemperatureScheduleGenerator():
	utility::VirtualBase(),
	max_temperature_(50.0),
	min_temperature_(0.62),
	current_round_(0),
	max_rounds_(0)
{}

/// @brief Options constructor.
HBP_TemperatureScheduleGenerator::HBP_TemperatureScheduleGenerator( core::Real const & max_temp, core::Real const & min_temp ) :
	utility::VirtualBase(),
	max_temperature_(max_temp),
	min_temperature_(min_temp),
	current_round_(0),
	max_rounds_(0)
{
	runtime_assert_string_msg( max_temp > 0.0, "Error in HBP_TemperatureScheduleGenerator::HBP_TemperatureScheduleGenerator(): The maximum temperature must be positive." );
	runtime_assert_string_msg( min_temp > 0.0, "Error in HBP_TemperatureScheduleGenerator::HBP_TemperatureScheduleGenerator(): The minimum temperature must be positive." );
}

/// @brief Copy constructor.
HBP_TemperatureScheduleGenerator::HBP_TemperatureScheduleGenerator( HBP_TemperatureScheduleGenerator const &src ) :
	VirtualBase( src ),
	max_temperature_( src.max_temperature_ ),
	min_temperature_( src.min_temperature_),
	current_round_( src.current_round_ ),
	max_rounds_( src.max_rounds_ )
{}

/// @brief Destructor.
HBP_TemperatureScheduleGenerator::~HBP_TemperatureScheduleGenerator(){}

///////////////////////////////////////// PUBLIC MEMBER FUNCTIONS /////////////////////////////////////////

/// @brief Set the maximum temperature for the run.
void
HBP_TemperatureScheduleGenerator::set_max_temperature(
	core::Real const & temp_in
) {
	runtime_assert_string_msg( temp_in > 0.0, "Error in HBP_TemperatureScheduleGenerator::set_max_temperature(): The maximum temperature must be positive." );
	max_temperature_ = temp_in;
}

/// @brief Set the minimum temperature for the run.
void
HBP_TemperatureScheduleGenerator::set_min_temperature(
	core::Real const & temp_in
) {
	runtime_assert_string_msg( temp_in > 0.0, "Error in HBP_TemperatureScheduleGenerator::set_min_temperature(): The minimum temperature must be positive." );
	min_temperature_ = temp_in;
}

/// @brief Set the index of the current round of simulated annealing.
void
HBP_TemperatureScheduleGenerator::set_current_round(
	core::Size const setting
) {
	if ( max_rounds_ > 0 ) runtime_assert( setting <= max_rounds_ );
	runtime_assert( setting > 0 );
	current_round_ = setting;
}

/// @brief Set the total nunmber of rounds of simulated annealing.
void
HBP_TemperatureScheduleGenerator::set_max_rounds(
	core::Size const setting
) {
	if ( current_round_ > 0 ) runtime_assert( current_round_ <= setting );
	runtime_assert( setting > 0 );
	max_rounds_ = setting;
}

/// @brief Get the maximum temperature, given the round.
/// @details Round must be set, and greater than zero.  The base class version of this function just
/// returns max_temperature_, but this can be overridden if a derived class wants to do something more
/// complicated.
core::Real
HBP_TemperatureScheduleGenerator::get_max_temp_for_current_round() const {
	return max_temperature_;
}

/// @brief Get the minimum temperature, given the round.
/// @details Round must be set, and greater than zero.  The base class version of this function just
/// returns min_temperature_, but this can be overridden if a derived class wants to do something more
/// complicated.
core::Real
HBP_TemperatureScheduleGenerator::get_min_temp_for_current_round() const {
	return min_temperature_;
}

} //protocols
} //helical_bundle_predict






