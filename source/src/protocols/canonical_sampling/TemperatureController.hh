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
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>

// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Protocol Headers
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/random/WeightedSampler.hh>
#include <boost/noncopyable.hpp>

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
/// simulation.
/// 
/// Methods are also get and set generally useful parameters.  These methods 
/// include get_high_temp(), get_low_temp(), get_temp_levels(), get_stride(), 
/// and temperature_level().  Temperature levels are provided for managing a 
/// number of discrete temperature levels.  Not all of these parameters are 
/// respected by and/or make sense for all temperature controllers.  For 
/// example, FixedTemperatureController ignores all of this information.

class TemperatureController : public protocols::moves::Mover {

// Public Interfaces
public:

	/// @brief Default constructor.
	TemperatureController();

	/// @brief Default constructor.
	TemperatureController(TemperatureController const & other);

	/// @brief No-op implemented only to satisfy the Mover interface.
	virtual
	void apply( core::pose::Pose& ) {}

	/// @brief Return false.  This class does not need to be reinitialized for 
	/// each job.
	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	/// @brief Return false.  This class does not need to be reinitialized for 
	/// new input.  
	virtual
	bool
	reinitialize_for_new_input() const { return false; }

	/// @brief Register the options used by this mover with the global options 
	/// system.
	static void register_options();

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	/// @copydoc ThermodynamicObserver::initialize_simulation
	virtual void
	initialize_simulation(
		core::pose::Pose & pose,
		MetropolisHastingsMover const & mover,
		core::Size cycle
	);

	/// @brief Execute a temperature move if necessary.
	/// @details This method is expected to return the new temperature (in units 
	/// of kT, although kT doesn't mean much in the context of rosetta).
	virtual
	core::Real
	temperature_move(
		core::pose::Pose & pose,
		MetropolisHastingsMover & mover,
		core::Real score
	)=0;

	/// @copydoc ThermodynamicObserver::finalize_simulation
	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		MetropolisHastingsMover const & mover
	);

// Getters and Setters
public:

	/// @brief Return the temperature of the underlying MonteCarlo object.
	core::Real temperature() const;

	/// @brief  Return the temperature of the given level.
	core::Real temperature( core::Size level ) const;

	/// @brief Return the current temperature level.
	/// @details Tempering controllers often work with a handful of discrete 
	/// temperature levels.  This method makes it possible to work with levels, 
	/// which are discrete, rather than temperatures, which are continuous.  
	/// Higher temperature levels correspond to higher temperatures, and as usual 
	/// in rosetta, counting starts at one.  So the lowest temperature being 
	/// simulated will always be level 1.
	/// @see n_temp_levels()
	/// @see temperature()
	core::Size temperature_level() const;

	/// @brief Return the number of temperature levels used by this controller.
	/// @details This parameter must be set from the command line using the 
	/// <tt>-tempering:temp:levels</tt> flag.  The default value is 10.  In some 
	/// cases (e.g. when using MpiParallelTempering), this parameter must also match 
	/// the number of processes allocated to your app by @c mpirun.  For example:
	/// @code{.sh}
	/// mpirun -np 32 my_parallel_tempering_app -tempering:temp:levels 32
	/// @endcode
	/// I'm not sure why this requirement exists. It's simple enough to ask MPI 
	/// how many threads are available; why not just set that many temperature 
	/// levels?  Perhaps this would make things harder to reproduce, because the 
	/// behavior of the program would depend on options passed to @c mpirun.
	/// @see temperature_level()
	/// @see temperature()
	core::Size n_temp_levels() const;

	/// @brief Return const access to the MonteCarlo object being controlled.
	protocols::moves::MonteCarloCOP
	monte_carlo() const;

	/// @brief Set the MonteCarlo object to be controlled.
	void
	set_monte_carlo(protocols::moves::MonteCarloOP monte_carlo);

// Protected Helpers
protected:

	/// @brief Return non-const access to the MonteCarlo object being controlled.
	protocols::moves::MonteCarloOP
	monte_carlo();

	/// @brief Help the constructor initialize new objects.
	void set_defaults();

	/// @brief Assign user-specified command-line values to data members.
	virtual
	void init_from_options();

	/// @brief Initialize temperatures and weights from a file.
	/// @details Return false if an IO error occurs.
	virtual
	bool init_from_file( std::string const& filename );

	/// @brief Save temperatures and weights to a file.
	virtual
	void write_to_file( std::string const& file_in, std::string const& output_name, utility::vector1< core::Real > const& wcounts );

	/// @brief Explicitly set the temperature levels by interpolating the given 
	/// parameters.
	void generate_temp_range( core::Real temp_low, core::Real temp_high, core::Size n_levels, InterpolationType interpolation = linear );

	/// @brief Explicitly set the temperature levels.
	void set_temperatures( utility::vector1< core::Real > const& );

	/// @brief Return the current temperature level.  Identical to 
	/// temperature_level() as far as I can tell.
	core::Size current_temp() { return current_temp_; }

	/// @brief Set the temperature to the given level.
	/// @details Note that the argument is a temperature level, not a raw 
	/// temperature.
	void set_current_temp( core::Size new_temp );

	/// @brief Assert that the current temperature of the MonteCarlo object 
	/// agrees with the current temperature level of this object.
	bool check_temp_consistency();

	/// @brief Forget all temperature levels and return to an uninitialized 
	/// state.
	void clear();

	/// @brief Return true if a temperature move should be made on this 
	/// iteration.
	bool time_for_temp_move() {
		return ++temp_trial_count_ % temperature_stride_ == 0;
	}

	/// @brief Return true if a statistics summary should be written.
	/// @see stats_silent_output()
	/// @see stats_file()
 	bool stats_line_output() const {
		return stats_line_output_;
	}

	/// @brief Return true if a statistics summary should be inserted into a 
	/// silent file.
	/// @see stats_line_output()
	/// @see stats_file()
 	bool stats_silent_output() const {
		return stats_silent_output_;
	}
	
	/// @brief Return the name of the silent file into which statistics should be 
	/// recorded.
	/// @see stats_line_output()
	/// @see stats_silent_output()
	std::string const& stats_file() const {
		return stats_file_;
	}

// Data Members
private:
	
	/// @brief The underlying monte carlo simulation.
	protocols::moves::MonteCarloOP monte_carlo_;

	/// @brief If false, the options used by this class have not yet been added 
	/// to the global options system.
	static bool options_registered_;

	/// @brief Temperature levels.
	utility::vector1< core::Real > temperatures_;

	/// @brief Frequency for attempting temperature moves.
	core::Size temperature_stride_;

	/// @brief If false, look for current temperature in monte_carlo_ before each 
	/// move.  Set to true by default.
	bool trust_current_temp_;

	/// @brief If true, a statistics summary will be written.
	bool stats_line_output_;

	/// @brief If true, the statistics summary will be inserted in a silent file.
	bool stats_silent_output_;

	/// @brief Name of the silent file used for writing statistics.
	std::string stats_file_;

	/// @brief Job object to report on temperatures.
	protocols::jd2::JobOP job_;

	/// @brief If false, init_from_options() will be called before the simulation 
	/// starts.
	bool instance_initialized_;

	/// @brief Current temperature level.  Not the current temperature!
	core::Size current_temp_;

	/// @brief Number of times time_for_temp_move() has been called.  This method 
	/// is meant to be called every time temperature_move() is called. 
	core::Size temp_trial_count_;

}; // TemperatureController

} // canonical_sampling
} // protocols

#endif
