// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBP_TemperatureScheduleGenerator.hh
/// @brief A base class to generate a temperature ramping scheme for a simulated annealing trajectory.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_helical_bundle_predict_HBP_TemperatureScheduleGenerator_hh
#define INCLUDED_protocols_helical_bundle_predict_HBP_TemperatureScheduleGenerator_hh

// Project headrs
#include <protocols/helical_bundle_predict/HBP_TemperatureScheduleGenerator.fwd.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace protocols {
namespace helical_bundle_predict {

/// @brief A base class to generate a temperature ramping scheme for a simulated annealing trajectory.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class HBP_TemperatureScheduleGenerator : public utility::VirtualBase {

public:

	/// @brief Default constructor.
	HBP_TemperatureScheduleGenerator();

	/// @brief Options constructor.
	HBP_TemperatureScheduleGenerator( core::Real const & max_temp, core::Real const & min_temp );

	/// @brief Copy constructor.
	HBP_TemperatureScheduleGenerator(HBP_TemperatureScheduleGenerator const & src);

	/// @brief Destructor.
	~HBP_TemperatureScheduleGenerator() override;

	/// @brief Copy this object and return an owning pointe to the copy.
	virtual HBP_TemperatureScheduleGeneratorOP clone() const = 0;

public: //Functions

	/// @brief Set the maximum temperature for the run.
	void set_max_temperature( core::Real const & temp_in );

	/// @brief Set the minimum temperature for the run.
	void set_min_temperature( core::Real const & temp_in );

	/// @brief Get the maximum temperature.
	inline core::Real const & max_temperature() const { return max_temperature_; }

	/// @brief Get the maximum temperature.
	inline core::Real const & min_temperature() const { return min_temperature_; }

	/// @brief Set the index of the current round of simulated annealing.
	void set_current_round( core::Size const setting );

	/// @brief Set the total nunmber of rounds of simulated annealing.
	void set_max_rounds( core::Size const setting );

	/// @brief Get the index of the current round of simulated annealing.
	inline core::Size current_round() const { return current_round_; }

	/// @brief Get the total nunmber of rounds of simulated annealing.
	inline core::Size max_rounds() const { return max_rounds_; }

	/// @brief Get the maximum temperature, given the round.
	/// @details Round must be set, and greater than zero.  The base class version of this function just
	/// returns max_temperature_, but this can be overridden if a derived class wants to do something more
	/// complicated.
	virtual core::Real get_max_temp_for_current_round() const;

	/// @brief Get the minimum temperature, given the round.
	/// @details Round must be set, and greater than zero.  The base class version of this function just
	/// returns min_temperature_, but this can be overridden if a derived class wants to do something more
	/// complicated.
	virtual core::Real get_min_temp_for_current_round() const;

	/// @brief Calculate the current temperature for this point in the trajectory.
	/// @details Pure virtual.  Must be implemented by derived classes.
	virtual core::Real calculate_current_temperature( core::Size const current_step, core::Size const max_steps ) const = 0;

private:

	/// @brief The maximum temperature for the run.
	/// @details Defaults to 50.
	core::Real max_temperature_;

	/// @brief The minimum temperature for the run.
	/// @details Defaults to 0.62.
	core::Real min_temperature_;

	/// @brief The index of the current round of simulated annealing.
	core::Size current_round_;

	/// @brief The total nunmber of rounds of simulated annealing.
	core::Size max_rounds_;

};


} //protocols
} //helical_bundle_predict



#endif //INCLUDED_protocols_helical_bundle_predict_HBP_TemperatureScheduleGenerator_hh





