// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file EnvClaim
/// @brief
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_claims_EnvClaim_hh
#define INCLUDED_protocols_environment_claims_EnvClaim_hh

// Unit Headers
#include <protocols/environment/claims/EnvClaim.fwd.hh>

// Package Headers
#include <core/environment/LocalPosition.fwd.hh>
#include <core/environment/FoldTreeSketch.hh>

#include <protocols/environment/claims/BrokerElements.hh>
#include <protocols/environment/ClaimingMover.fwd.hh>
#include <protocols/environment/ProtectedConformation.fwd.hh>

// Project Headers
#include <core/types.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>
#include <utility/vector1.hh>


// option key includes


namespace protocols {
namespace environment {
namespace claims {
/// A better EnvClaims class would provide some extracting functions:
/// by owner
/// by type

class EnvClaim : public utility::pointer::ReferenceCount {

  typedef core::environment::FoldTreeSketch FoldTreeSketch;

public:
  ///@brief Virtual destructor
  virtual ~EnvClaim();

  EnvClaim( ClaimingMoverOP );

  virtual EnvClaimOP clone() const = 0;

  ClaimingMoverOP owner() const;

  ///@brief build ResidueElements that indicate the introduction of a new peptide edge into the fold tree.
  virtual void yield_elements( FoldTreeSketch const&, ResidueElements& ) const {};

  ///@brief build the JumpElements that represent the inclusion of a jump in the nascent FoldTree
  virtual void yield_elements( FoldTreeSketch const&, JumpElements& ) const {};

  ///@brief build and export the CutElements that represent the inclusion of a cut in the tree.
  virtual void yield_elements( FoldTreeSketch const&, CutElements& ) const {};

  ///@brief build and export the CutElements that represent the inclusion of a cut in the tree.
  virtual void yield_elements( FoldTreeSketch const&, CutBiasElements& ) const {};

  ///@brief build and export DOFElements, which represent control over non-jump dofs (torsions, bond lengths, angles) final conformation.
  virtual void yield_elements( ProtectedConformationCOP const&, DOFElements& ) const {};

  virtual std::string str_type() const = 0;

  virtual void show( std::ostream& os ) const;

protected:
  virtual
  DOFElement wrap_dof_id( core::id::DOF_ID const& id ) const;

private:

  ClaimingMoverOP claim_source_;

}; //class EnvClaim

extern std::ostream& operator<<( std::ostream& os, EnvClaim const& );
extern std::ostream& operator<<( std::ostream& os, EnvClaims const& );

}
}
}

#endif
