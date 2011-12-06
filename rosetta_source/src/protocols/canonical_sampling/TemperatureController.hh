// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/canonical_sampling/MetropolisHastingsMover.hh
/// @brief
/// @author Oliver Lange ( oliver.lange@tum.de )

#ifndef INCLUDED_protocols_canonical_sampling_TemperatureController_hh
#define INCLUDED_protocols_canonical_sampling_TemperatureController_hh

// Unit Headers
#include <protocols/canonical_sampling/TemperatureController.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/WeightedSampler.hh>
#include <protocols/jd2/Job.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace canonical_sampling {

typedef enum { linear, exponential } InterpolationType;

std::string
interpolation_type_enum_to_string( InterpolationType interp_enum );

InterpolationType
interpolation_type_string_to_enum( std::string const & interp_string );

///@details
class TemperatureController : public protocols::canonical_sampling::ThermodynamicObserver {

public:

	TemperatureController();

	TemperatureController( TemperatureController const& );

	virtual
	void apply( core::pose::Pose& ) {};

	virtual
	std::string
	get_name() const;

	virtual
	bool
	reinitialize_for_each_job() const { return false; };

	virtual
	bool
	reinitialize_for_new_input() const { return false; };

	virtual
	void
	observe_after_metropolis(
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	);
	/// @brief execute the temperatur move ( called by observer_after_metropolis )
	/// returns the current temperatur in kT.
	virtual core::Real temperature_move( core::Real score) = 0;

	virtual core::Real temperature() const = 0;

	///@brief  return temperature of a certain level
	virtual core::Real temperature( core::Size level ) const = 0 ;


	virtual core::Size temperature_level() const {
		return 1;
	}

	virtual core::Size n_temp_levels() const { return 1; };

	protocols::moves::MonteCarloCOP
	monte_carlo() const;

	void
	set_monte_carlo(
		protocols::moves::MonteCarloOP monte_carlo
	);

protected:
	protocols::moves::MonteCarloOP
	monte_carlo();

private:
	protocols::moves::MonteCarloOP monte_carlo_;
}; //end TemperatureController


class FixedTemperatureController : public protocols::canonical_sampling::TemperatureController {
public:
	FixedTemperatureController( core::Real temp ) :
		temperature_ ( temp ) {};

	protocols::moves::MoverOP
	clone() const { return new protocols::canonical_sampling::FixedTemperatureController( temperature_ ); };

	virtual
	std::string
	get_name() const { return "FixedTemperatureContoller"; }
	/// @brief execute the temperatur move ( called by observer_after_metropolis )
	/// returns the current temperatur in kT.
	virtual core::Real temperature_move( core::Real score) { return temperature_; };

	virtual core::Real temperature() const { return temperature_; };

	///@brief  return temperature of a certain level
	virtual core::Real temperature( core::Size level ) const { return temperature_; };

private:
	core::Real temperature_;
};

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_TemperatureController_HH
