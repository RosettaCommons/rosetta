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


#ifndef INCLUDED_protocols_environment_claims_CutBiasClaim_hh
#define INCLUDED_protocols_environment_claims_CutBiasClaim_hh


// Unit Headers
#include <protocols/environment/claims/CutBiasClaim.fwd.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>
#include <core/environment/FoldTreeSketch.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/BrokerElements.hh>

// Project Headers
#include <core/id/types.hh>
#include <core/fragment/SecondaryStructure.hh>

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

class CutBiasClaim : public EnvClaim {
  typedef EnvClaim Parent;
  typedef core::environment::FoldTreeSketch FoldTreeSketch;

public:
  // Initializer for a single backbone angle
  CutBiasClaim( ClaimingMoverOP owner,
                std::string const& label,
                core::fragment::SecondaryStructure const& ss );

  virtual void yield_elements( FoldTreeSketch const& fts, CutBiasElements& elements ) const;

  EnvClaimOP clone() const;

  virtual std::string str_type() const;

  virtual void show( std::ostream& os ) const;

private:
  std::string label_;
  core::fragment::SecondaryStructure ss_;


}; //class CutBiasClaim


}
}
}

#endif
