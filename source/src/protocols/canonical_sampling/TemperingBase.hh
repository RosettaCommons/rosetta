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
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/WeightedSampler.hh>
#include <protocols/jd2/Job.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace canonical_sampling {

///@details
class TemperingBase : public protocols::canonical_sampling::TemperatureController {
	typedef TemperatureController Parent;
public:

	TemperingBase();

	TemperingBase( TemperingBase const& );

	virtual
	void apply( core::pose::Pose& ) {};

	virtual
	std::string
	get_name() const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	/// @brief callback executed before any Monte Carlo trials
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

	/// @brief callback executed after all Monte Carlo trials
	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	/// @brief return current_temperature (in monte-carlo object)
	core::Real temperature() const;

	///@brief  return temperature of a certain level
	core::Real temperature( core::Size level ) const;

	core::Size temperature_level() const {
		return current_temp_;
	}

	/// @brief callback executed after all Monte Carlo trials
	core::Size n_temp_levels() const;

protected:
	void set_defaults();

	/// @brief Assigns user specified values to primitive members using command line options
	virtual
	void init_from_options();

	/// @brief initialize temperatures and weights from file, return false if IO error occurrs
	virtual
	bool initialize_from_file( std::string const& filename );

	virtual
	void write_to_file( std::string const& file_in, std::string const& output_name, utility::vector1< core::Real > const& wcounts );

	bool check_temp_consistency();

	bool time_for_temp_move() {
		return ++temp_trial_count_ % temperature_stride_ == 0;
	}

protected:

	core::Size current_temp() {
		return current_temp_;
	}

	//back to uninitialized state
	void clear();

	void set_temperatures( utility::vector1< core::Real > const& );

	void set_current_temp( core::Size new_temp );

 	bool stats_line_output() const {
		return stats_line_output_;
	}
 	bool stats_silent_output() const {
		return stats_silent_output_;
	}
	std::string const& stats_file() const {
		return stats_file_;
	}

	void generate_temp_range( core::Real temp_low, core::Real temp_high, core::Size n_levels, InterpolationType interpolation = linear );

/// ------------------ register cmdline options ---------------------------
private:
	static bool options_registered_;

public:
	static void register_options();

/// ---------------- member variables --------------------------
private:
	///---  configurables... -----
	//temperature levels
	utility::vector1< core::Real > temperatures_;

	//attempt frequency for temperature moves
	core::Size temperature_stride_;

	//if false look for current temperature in monte_carlo_ at each move (default true)
	bool trust_current_temp_;

	//how should statistics output be written -- common options read by child-classes
	bool stats_line_output_;
	bool stats_silent_output_;
	std::string stats_file_;

	//job object to report on temperatures
	protocols::jd2::JobOP job_;

	//if not initialized when simulations starts call init_from_options()
	bool instance_initialized_;

	//current temperature level
	core::Size current_temp_;

	//counting calls to temp_moves (-->temperature_stride_)
	core::Size temp_trial_count_;
}; //end TemperingBase

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_TemperingBase_HH
