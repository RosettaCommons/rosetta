// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange


#ifndef INCLUDED_protocols_environment_claims_TorsionClaim_hh
#define INCLUDED_protocols_environment_claims_TorsionClaim_hh


// Unit Headers
#include <protocols/environment/claims/TorsionClaim.fwd.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>
#include <core/environment/FoldTreeSketch.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/BrokerElements.hh>

// Project Headers
#include <core/id/types.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//#include <basic/options/option_macros.hh>

//// C++ headers
#include <string>
#include <sstream>
#include <utility>

namespace protocols {
namespace environment {
namespace claims {

class TorsionClaim : public EnvClaim {
  typedef core::environment::FoldTreeSketch FoldTreeSketch;

public:
  typedef core::environment::LocalPosition LocalPosition;
  typedef core::environment::LocalPositions LocalPositions;

  // Initializer for a single backbone angle
  TorsionClaim( ClaimingMoverOP owner,
                LocalPosition const& local_pos );

  // Initializer for a contiguous range of backbone angles.
  TorsionClaim( ClaimingMoverOP owner,
                std::string const& label,
                std::pair< core::Size, core::Size > const& range );

  virtual void yield_elements( FoldTreeSketch const& fts, TorsionElements& elements ) const;

  LocalPositions const& positions() const;

  ControlStrength const& ctrl_strength() const;

  void ctrl_strength( ControlStrength const& );

  InitializationStrength const& init_strength() const;

  void init_strength( InitializationStrength const& );

  virtual EnvClaimOP clone() const;

  virtual std::string str_type() const;

  virtual void show( std::ostream& os ) const;

private:
  LocalPositions local_positions_;
  ControlStrength c_str_;
  InitializationStrength i_str_;

}; //class TorsionClaim


}
}
}

#endif
