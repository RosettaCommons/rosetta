// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/canonical_sampling/MultiTemperatureTrialCounter.hh
/// @brief
/// @author Oliver Lange

#ifndef INCLUDED_protocols_canonical_sampling_MultiTemperatureTrialCounter_hh
#define INCLUDED_protocols_canonical_sampling_MultiTemperatureTrialCounter_hh

// Unit Headers
#include <protocols/canonical_sampling/TemperatureController.hh>
#include <protocols/moves/TrialCounter.hh>

// Project Headers
// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace canonical_sampling {

class MultiTemperatureTrialCounter {
public:
  MultiTemperatureTrialCounter() {};
  MultiTemperatureTrialCounter( TemperatureControllerCOP );

  void reset();

  void count_trial( std::string const& );
  void count_accepted( std::string const& );
  void count_energy_drop( std::string const&, core::Real );

  // Undefined, commenting out to fix PyRosetta build  core::Size trial( std::string const& );
  // Undefined, commenting out to fix PyRosetta build  core::Size accepted( std::string const& );
  // Undefined, commenting out to fix PyRosetta build  core::Real energy_drop( std::string const& );

  protocols::moves::TrialCounter const& operator[]( core::Size ) const;
  protocols::moves::TrialCounter& operator[]( core::Size );

  void show( std::ostream& ) const;
  void show() const;
  void write_to_file( std::string const& file, std::string const& tag ) const;

  void set_temperature_observer( TemperatureControllerCOP );

private:
	void _write_to_stream( std::ostream&, std::string const& tag ) const;
  TemperatureControllerCOP tempering_;
  utility::vector1< protocols::moves::TrialCounter > counters_;
};


}
}

#endif
