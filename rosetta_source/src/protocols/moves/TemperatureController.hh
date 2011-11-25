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

/// @file /protocols/moves/MetropolisHastingsMover.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_moves_TemperatureController_hh
#define INCLUDED_protocols_moves_TemperatureController_hh

// Unit Headers
#include <protocols/moves/TemperatureController.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/ThermodynamicMover.hh>
#include <protocols/moves/ThermodynamicObserver.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/WeightedSampler.hh>
#include <protocols/jd2/Job.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace moves {

///@details
class TemperatureController : public ThermodynamicObserver {

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

	/// @brief callback executed before any Monte Carlo trials
	virtual
	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	virtual
	void
	observe_after_metropolis(
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
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

	virtual/// @brief callback executed before any Monte Carlo trials
	void
	initialize_simulation() {};

	virtual core::Size n_temp_levels() const { return 1; };

	MonteCarloCOP
	monte_carlo() const;

	void
	set_monte_carlo(
		MonteCarloOP monte_carlo
	);

private:
	MonteCarloOP monte_carlo_;
}; //end TemperatureController


class FixedTemperatureController : public TemperatureController {
public:
	FixedTemperatureController( core::Real temp ) :
		temperature_ ( temp ) {};

	MoverOP
	clone() const { new FixedTemperatureController( temperature_ ); };

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

} //namespace moves
} //namespace protocols

#endif //INCLUDED_protocols_moves_TemperatureController_HH
