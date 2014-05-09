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

/// @brief Base class for controlling the temperature of a simulation.
///
/// @details Many schemes to improve the performance of condensed phase Monte 
/// Carlo simulations depends on changing the temperature of the system.  
/// Common examples include simulated annealing and parallel tempering.  This 
/// class provides an interface for writing these algorithms.  The most 
/// important method is temperature_move(), which is responsible for actually 
/// changing the temperature of the MonteCarlo object used for the underlying 
/// simulation.  Methods like temperature_level() are also provided for 
/// managing a discrete number of different temperature levels, which is a 
/// common feature of these algorithms.
/// 
/// The TemperingBase class serves a similar role to this one, but is geared 
/// towards controllers that actually intend to change the temperature.  This 
/// class also is parent to FixedTemperatureController, which is the default 
/// controller used by MetropolisHastingsMover.

class TemperatureController : public ThermodynamicObserver {

public:

	/// @brief Default constructor.
	TemperatureController();

	/// @brief Copy constructor.
	TemperatureController( TemperatureController const& );

	/// @brief No-op implemented only to satisfy the Mover interface.
	virtual
	void apply( core::pose::Pose& ) {};

	/// @brief Return the name of this class.
	virtual
	std::string
	get_name() const;

	/// @brief Return false.  This class does not need to be reinitialized for 
	/// each job.
	virtual
	bool
	reinitialize_for_each_job() const { return false; };

	/// @brief Return false.  This class does not need to be reinitialized for 
	/// new input.  
	virtual
	bool
	reinitialize_for_new_input() const { return false; };

	virtual
	void
	observe_after_metropolis(
		MetropolisHastingsMover const & metropolis_hastings_mover
	);

	/// @brief Execute the temperature move.
	/// @details This method is called by observe_after_metropolis() and is 
	/// expected to return the new temperature (in units of kT, to the extent 
	/// that that is meaningful in the context of rosetta).
	virtual
	core::Real
	temperature_move( core::Real score ) = 0;

	/// @brief Execute a temperature move which depends on the current pose.
	/// @details The default implementation just calls the pose-independent 
	/// temperature_pose() method with the energy of the given pose.  However, 
	/// the HamiltonianExchange temperature controller needs to evaluate the 
	/// alternative Hamiltonian.
	virtual
	core::Real
	temperature_move( core::pose::Pose& pose );

	/// @brief Return the current temperature.
	virtual core::Real temperature() const = 0;

	/// @brief Set the current temperature to match given level.
	virtual core::Real temperature( core::Size level ) const = 0 ;

	/// @brief Return the current temperature level.
	/// @details Tempering controllers often work with a handful of discrete 
	/// temperature levels.  This method makes it possible to work with levels, 
	/// which are discrete, rather than temperatures, which are continuous.
	/// @see n_temp_levels()
	/// @see temperature()
	virtual core::Size temperature_level() const {
		return 1;
	}

	virtual
	void
	initialize_simulation(
		core::pose::Pose & pose,
		MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size level,
		core::Real temperature,
		core::Size cycle
	);

	virtual
	void
	initialize_simulation(
		core::pose::Pose &,
		MetropolisHastingsMover const &,
		core::Size
	) {};

	/// @brief Return the number of temperature levels used by this controller.
	/// @see temperature_level()
	virtual core::Size n_temp_levels() const { return 1; };

	/// @brief Return const access to the MonteCarlo object being controlled.
	protocols::moves::MonteCarloCOP
	monte_carlo() const;

	/// @brief Return true if the simulation has been completed.
	virtual
	bool finished_simulation( core::Size trials, core::Size ntrials) {
		return trials > ntrials;
	}

	/// @brief Set the MonteCarlo object to be controlled.
	virtual
	void set_monte_carlo(
		protocols::moves::MonteCarloOP monte_carlo
	);

protected:
	/// @brief Return non-const access to the MonteCarlo object being controlled.
	protocols::moves::MonteCarloOP
	monte_carlo();

private:
	protocols::moves::MonteCarloOP monte_carlo_;
}; //end TemperatureController


/// @brief Maintain a constant temperature.
/// @details This is the default temperature controller used by 
/// MetropolisHastingsMover.
class FixedTemperatureController : public TemperatureController {
public:

	/// @brief Constructor with temperature parameter.
	FixedTemperatureController( core::Real temp ) :
		temperature_ ( temp ) {};

	/// @brief Return a copy of this mover.
	protocols::moves::MoverOP
	clone() const { return new protocols::canonical_sampling::FixedTemperatureController( temperature_ ); };

	virtual
	std::string
	get_name() const { return "FixedTemperatureContoller"; }

	/// @brief Return the same constant temperature every time.
	virtual core::Real temperature_move( core::Real ) { return temperature_; };

	virtual core::Real temperature() const { return temperature_; };

	virtual core::Real temperature( core::Size ) const { return temperature_; };

private:
	core::Real temperature_;
};

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_TemperatureController_HH
