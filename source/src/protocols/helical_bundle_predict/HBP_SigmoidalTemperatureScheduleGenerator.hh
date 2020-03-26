// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBP_SigmoidalTemperatureScheduleGenerator.hh
/// @brief A class to generate a temperature ramping scheme for a simulated annealing trajectory.  This version
/// ramps sigmoidally.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_helical_bundle_predict_HBP_SigmoidalTemperatureScheduleGenerator_hh
#define INCLUDED_protocols_helical_bundle_predict_HBP_SigmoidalTemperatureScheduleGenerator_hh

// Project headrs
#include <protocols/helical_bundle_predict/HBP_SigmoidalTemperatureScheduleGenerator.fwd.hh>
#include <protocols/helical_bundle_predict/HBP_TemperatureScheduleGenerator.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace protocols {
namespace helical_bundle_predict {

/// @brief A class to generate a temperature ramping scheme for a simulated annealing trajectory.  This version
/// ramps sigmoidally.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class HBP_SigmoidalTemperatureScheduleGenerator : public HBP_TemperatureScheduleGenerator {

public:

	/// @brief Default constructor.
	HBP_SigmoidalTemperatureScheduleGenerator();

	/// @brief Options constructor.
	HBP_SigmoidalTemperatureScheduleGenerator( core::Real const & max_temp, core::Real const & min_temp );

	/// @brief Copy constructor.
	HBP_SigmoidalTemperatureScheduleGenerator(HBP_SigmoidalTemperatureScheduleGenerator const & src);

	/// @brief Destructor.
	~HBP_SigmoidalTemperatureScheduleGenerator() override;

	/// @brief Copy this object and return an owning pointe to the copy.
	HBP_TemperatureScheduleGeneratorOP clone() const override;

public: //Functions

	/// @brief Calculate the current temperature for this point in the trajectory.
	core::Real calculate_current_temperature( core::Size const current_step, core::Size const max_steps ) const override;

	/// @brief Set the exponent.
	/// @details This is a value that determines how S-shaped the sigmoid is.  A value of 1 gives a
	/// straight line; higher values have more lag at start and end.
	/// @note Must be greater than zero.
	void set_exponent( core::Real const & setting );

private:

	/// @brief Value that determines how S-shaped the sigmoid is.  A value of 1 gives a
	/// straight line; higher values have more lag at start and end.
	/// @details Defaults to 2.0.
	core::Real exponent_;

};


} //protocols
} //helical_bundle_predict



#endif //INCLUDED_protocols_helical_bundle_predict_HBP_SigmoidalTemperatureScheduleGenerator_hh





