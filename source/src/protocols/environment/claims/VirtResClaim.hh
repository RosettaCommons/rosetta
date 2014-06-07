// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/claims/VirtResClaim.hh
/// @brief A claim for expressing the desire for new sequence to be inserted at a particular point.
/// @author Justin Porter


#ifndef INCLUDED_protocols_environment_claims_VirtResClaim_hh
#define INCLUDED_protocols_environment_claims_VirtResClaim_hh

// Unit Headers
#include <protocols/environment/claims/VirtResClaim.fwd.hh>
#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/JumpClaim.hh>
#include <protocols/environment/claims/XYZClaim.hh>

// Package Headers
#include <core/environment/FoldTreeSketch.hh>
#include <core/environment/LocalPosition.hh>

#include <protocols/environment/claims/BrokerElements.hh>

// Project Headers

// ObjexxFCL Headers

// Utility headers

//// C++ headers

// option key includes


namespace protocols {
namespace environment {
namespace claims {

class VirtResClaim : public EnvClaim {
  typedef EnvClaim Parent;
  typedef core::environment::LocalPosition LocalPosition;
  typedef core::environment::FoldTreeSketch FoldTreeSketch;

public:

  VirtResClaim( ClaimingMoverOP owner,
                LocalPosition parent,
                std::string const& jump_label,
                std::string const& vrt_label );

  virtual void yield_elements( FoldTreeSketch const& fts, ResidueElements& elements ) const;

  virtual void yield_elements( FoldTreeSketch const& fts, JumpElements& elements ) const;

  virtual void yield_elements( FoldTreeSketch const& fts, CutElements& elements ) const;

  virtual void yield_elements( ProtectedConformationCOP const&, DOFElements& elements ) const;

  JumpClaim& jump() { return j_claim_; }

  void strength( ControlStrength const& cstr, ControlStrength const& istr );

  std::string const& vrt_label() const;

  virtual EnvClaimOP clone() const;

  virtual std::string str_type() const;

  virtual void show( std::ostream& os ) const;

private:
  std::string vrt_label_;
  JumpClaim j_claim_;
  XYZClaim xyz_claim_;

}; //VirtResClaim

}
}
}

#endif
