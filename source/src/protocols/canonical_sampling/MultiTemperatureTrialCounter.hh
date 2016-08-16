// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/canonical_sampling/MultiTemperatureTrialCounter.hh
/// @brief
/// @author Oliver Lange

#ifndef INCLUDED_protocols_canonical_sampling_MultiTemperatureTrialCounter_hh
#define INCLUDED_protocols_canonical_sampling_MultiTemperatureTrialCounter_hh

// Unit Headers
#include <protocols/canonical_sampling/TemperatureController.hh>
#include <protocols/moves/TrialCounter.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace canonical_sampling {

/// @brief Keep track of trial statistics for any number of replicas.
/// @details This class helps MetropolisHastingsMover keep track of move
/// statistics.  At the end of a simulation, operator[]() can be used to access
/// the TrialCounter objects kept for each temperature level.  Alternatively,
/// the show() and write_to_file() methods can also be used to directly output
/// acceptance rates to a stream or file.
class MultiTemperatureTrialCounter {
public:

	/// @brief Default constructor.  A temperature controller must be set before
	/// the trial counter can be used.
	MultiTemperatureTrialCounter() {};

	/// @brief Fully construct the counter with a temperature controller.
	MultiTemperatureTrialCounter( TemperatureController const * );

	/// @brief Set all counters for all temperatures to zero.
	void reset();

	/// @brief Note that a move of the given type was attempted.
	void count_trial( std::string const& );

	/// @brief Note that a move of the given type was accepted.
	void count_accepted( std::string const& );

	/// @brief Note that a move of the given type led to the given energy drop.
	void count_energy_drop( std::string const&, core::Real );

	/// @brief Return const access to the TrialCounter for the given temperature
	/// level.
	protocols::moves::TrialCounter const& operator[]( core::Size ) const;

	/// @brief Return non-const access to the TrialCounter for the given
	/// temperature level.
	protocols::moves::TrialCounter& operator[]( core::Size );

	/// @brief Write acceptance rates for each move at each temperature to the
	/// given stream.
	void show( std::ostream& ) const;

	/// @brief Write acceptance rates for each move at each temperature to this
	/// module's tracer.
	void show() const;

	/// @brief Write acceptance rates for each move at each temperature to the
	/// given file.
	void write_to_file( std::string const& file, std::string const& tag ) const;

	/// @brief Set the temperature controller.
	void set_temperature_observer( TemperatureController const * );

private:

	/// @brief Help write the trial counters to the given stream.
	/// @see show()
	/// @see write_to_file()
	void _write_to_stream( std::ostream&, std::string const& tag ) const;

	TemperatureController const * tempering_;
	utility::vector1< protocols::moves::TrialCounter > counters_;
};


}
}

#endif
