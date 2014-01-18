// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/moves/TrialCounter.hh
/// @brief
/// @author Oliver Lange

#ifndef INCLUDED_protocols_moves_TrialCounter_hh
#define INCLUDED_protocols_moves_TrialCounter_hh

// Unit Headers
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/WeightedSampler.hh>
#include <protocols/jd2/Job.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>

namespace protocols {
namespace moves {

class TrialCounter {
public:
  TrialCounter() {};
  void reset();

  void count_trial( std::string const& );
  void count_accepted( std::string const& );
  void count_energy_drop( std::string const&, core::Real );

  core::Size trial( std::string const& ) const;
  core::Size accepted( std::string const& ) const;
  core::Real energy_drop( std::string const& ) const;
	utility::vector1< std::string > const tags () const;

  void show( std::ostream&, std::string line_header="", bool with_end_line = true ) const;
  void show() const;
  core::Size total_trials() const;
private:
  std::map< std::string, int > trial_counter_;
  std::map< std::string, int > accept_counter_;
  std::map< std::string, core::Real > energy_drop_counter_;
};


}
}

#endif
