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

// Unit headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/TrialCounter.fwd.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>

namespace protocols {
namespace moves {

class TrialCounter : public utility::pointer::ReferenceCount {
public:
  TrialCounter() {};
  virtual void reset();

  virtual void count_trial( std::string const& );
  virtual void count_accepted( std::string const& );
  virtual void count_energy_drop( std::string const&, core::Real );

  virtual core::Size total_trials() const;
  virtual core::Size trial( std::string const& ) const;
  virtual core::Size accepted( std::string const& ) const;
  virtual core::Real energy_drop( std::string const& ) const;
	virtual utility::vector1< std::string > const tags () const;

  void show() const;
  virtual void show(
			std::ostream&,
			std::string line_header="",
			bool with_end_line = true ) const;

private:
  std::map< std::string, int > trial_counter_;
  std::map< std::string, int > accept_counter_;
  std::map< std::string, core::Real > energy_drop_counter_;
};


}
}

#endif
