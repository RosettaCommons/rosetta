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
/// @author  Oliver Lange ( oliver.lange@tum.de )

#ifndef INCLUDED_protocols_canonical_sampling_TemperingBase_hh
#define INCLUDED_protocols_canonical_sampling_TemperingBase_hh

// Unit Headers
#include <protocols/canonical_sampling/TemperingBase.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/canonical_sampling/TemperatureController.hh>
#include <protocols/canonical_sampling/MultiTemperatureTrialCounter.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/WeightedSampler.hh>
#include <protocols/jd2/Job.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace canonical_sampling {

/// @brief Base class for tempering Monte Carlo optimizations.
///
/// @details Many important Monte Carlo optimization techniques, like simulated
/// annealing and parallel tempering, depend on a changing temperature
/// schedule.  TemperatureController provides the essential interface for
/// providing this functionality.  This class provides a lot of useful
/// protected member functions, especially with regard to input (i.e.
/// command-line or file) and output (i.e. tracer or silent file).
///
/// That said, my first impression is that this class really limits what you
/// can do in terms of output.  The emphasis seems to be on silent files, so it
/// would be hard to instead use (for example) database output.  In general,
/// putting IO code in a base class seems like a bad idea.  Better to do that
/// kind of stuff with object composition, so different IO formats can easily
/// be swapped in and out.  Perhaps this would be a good target for a small
/// refactoring project.

class TemperingBase : public protocols::canonical_sampling::TemperatureController {
	typedef TemperatureController Parent;
public:

	/// @brief Default constructor.
	TemperingBase();

	/// @brief Copy constructor.
	TemperingBase( TemperingBase const& );

	/// @brief No-op implemented only to satisfy the Mover interface.
	virtual
	void apply( core::pose::Pose& ) {};

	virtual
	std::string
	get_name() const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	virtual void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	);

	virtual
	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size level,
		core::Real temperature,
		core::Size cycle
	);

	virtual
	void
	observe_after_metropolis(
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	/// @brief Return the temperature of the underlying MonteCarlo object.
	core::Real temperature() const;

	/// @brief  Return the temperature of the given level.
	core::Real temperature( core::Size level ) const;

	core::Size temperature_level() const {
		return current_temp_;
	}

	core::Size n_temp_levels() const;

protected:
	/// @brief Help the constructor initialize the object.
	void set_defaults();

	/// @brief Assign user-specified command-line values to data members.
	virtual
	void init_from_options();

	/// @brief Initialize temperatures and weights from a file.
	/// @details Return false if an IO error occurs.
	virtual
	bool initialize_from_file( std::string const& filename );

	/// @brief Save temperatures and weights to a file.
	virtual
	void write_to_file( std::string const& file_in, std::string const& output_name, utility::vector1< core::Real > const& wcounts );

	/// @brief Assert that the current temperature of the MonteCarlo object
	/// agrees with the current temperature level of this object.
	bool check_temp_consistency();

	/// @brief Return true if a temperature move should be made on this
	/// iteration.
	virtual
	bool time_for_temp_move() {
		return ++temp_trial_count_ % temperature_stride_ == 0;
	}

	void reset_temp_counter() {
		temp_trial_count_ = 0;
	}

protected:

	/// @brief Return the current temperature level.  Identical to
	/// temperature_level() as far as I can tell.
	core::Size current_temp() const {
		return current_temp_;
	}

	/// @brief Forget all temperature levels and return to an uninitialized
	/// state.
	void clear();

	/// @brief Explicitly set the temperature levels.
	void set_temperatures( utility::vector1< core::Real > const& );

	/// @brief Set the temperature to the given level.
	/// @details Note that the argument is a temperature level, not a raw
	/// temperature.
	virtual
	void set_current_temp( core::Size new_temp );

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

	/// @brief Explicitly set the temperature levels by interpolating the given
	/// parameters.
	void generate_temp_range( core::Real temp_low, core::Real temp_high, core::Size n_levels, InterpolationType interpolation = linear );

	MultiTemperatureTrialCounter& trial_counter() {
		return trial_counter_;
	}

private:
	static bool options_registered_;

public:
	/// @brief Register the options used by this mover with the global options
	/// system.
	static void register_options();

private:
protected:

	/// @brief Temperature levels.
	utility::vector1< core::Real > temperatures_;

	/// @brief Frequency for attempting temperature moves (e.g. once every
	/// `io_stride_` steps).
	core::Size temperature_stride_;

	/// @brief Frequency with which statistics should be written (e.g. once every
	/// `io_stride_` steps).
	core::Size io_stride_;

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

	MultiTemperatureTrialCounter trial_counter_;

}; //end TemperingBase

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_TemperingBase_HH
