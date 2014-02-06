// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/canonical_sampling/MultiTempTrialCounter.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_MultiTempTrialCounter_hh
#define INCLUDED_protocols_canonical_sampling_MultiTempTrialCounter_hh

// Unit headers
#include <protocols/canonical_sampling/MultiTempTrialCounter.fwd.hh>
#include <protocols/canonical_sampling/TemperatureController.fwd.hh>

// Protocol headers
#include <protocols/moves/TrialCounter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace canonical_sampling {
 
/// @brief Keep track of move statistics at different temperature levels.
/// 
/// @details Moves are more likely to be accepted at higher temperatures, so it 
/// doesn't make sense to lump all the moves made at different temperatures 
/// into the same statistics.  Since the default TrialCounter used by 
/// MonteCarlo is agnostic to the temperature, this class is provided as a 
/// drop-in replacement.  Use MonteCarlo.set_counter() to replace the default 
/// counter with a MultiTempTrialCounter.  Once the simulation is complete, use 
/// temp_level() to return the move statistics for the given temperature level.  
/// The trial(), accepted(), and energy_drop() method are provided to satisfy 
/// the TrialCounter interface, but are not meaningful because they merge the 
/// results from all the different temperatures.

class MultiTempTrialCounter : public protocols::moves::TrialCounter {

public:

	/// @brief Default constructor.
	/// @details If the number of temperature levels managed by the controller 
	/// changes after this object is constructed, reset() must be called.
	MultiTempTrialCounter(TemperatureControllerCOP);

	/// @brief Default destructor.
	~MultiTempTrialCounter();

	/// @brief Set all counters for all temperatures to zero.
  void reset();

	/// @brief Set the temperature controller.
  void reset_temp_controller(TemperatureControllerCOP);

	/// @brief Record when a move of the given type was attempted.
	void count_trial(std::string const & tag);

	/// @brief Record when a move of the given type was accepted.
	void count_accepted(std::string const & tag);

	/// @brief Record when a move of the given type led to the given energy drop.
	void count_energy_drop(std::string const & tag, core::Real drop);

	/// @brief Copy all statistics to the root node, if MPI is enabled.
	void collect();

	/// @brief Return const access to the TrialCounter for the given temperature 
	/// level.
  protocols::moves::TrialCounter const & temp_level(core::Size) const;

	/// @brief Return the number of temperature levels recorded by this observer.
	core::Size num_temp_levels() const;

	/// @brief Return how many trials have been recorded for all moves and all 
	/// temperatures.
	core::Size total_trials() const;

	/// @brief Return how many times the given move was attempted for all 
	/// temperatures.
	core::Size trial(std::string const & tag) const;

	/// @brief Return how many times the given move was accepted for all 
	/// temperatures.
	core::Size accepted(std::string const & tag) const;

	/// @brief Return the average drop in energy effected by this move for all 
	/// temperatures.
	core::Real energy_drop(std::string const & tag) const;

	/// @brief Return the tag name for any move recorded at any temperature.
	virtual utility::vector1< std::string > const tags() const;

	/// @brief Write acceptance rates for each move at each temperature to the 
	/// given stream.
	void show(
			std::ostream&,
			std::string line_header="",
			bool with_end_line=true) const;

	/// @brief Write acceptance rates for each move at each temperature to the 
	/// given file.
  void write_to_file( std::string const& file, std::string const& tag ) const;

private:

	/// @brief Return a reference to the TrialCounter representing the current 
	/// temperature level.  Some sanity checks are done in debug mode.
	/// @see count_trial()
	/// @see count_accepted()
	/// @see count_energy_drop()
  protocols::moves::TrialCounter & get_current_counter();

	/// @brief Help write the trial counters to the given stream.
	/// @see show()
	/// @see write_to_file()
	void write_to_stream( std::ostream&, std::string const& tag ) const;

private:
  TemperatureControllerCOP temp_controller_;
  utility::vector1< protocols::moves::TrialCounter > counters_;

}; // end MultiTempTrialCounter

} // end canonical_sampling
} // end protocols

#endif
