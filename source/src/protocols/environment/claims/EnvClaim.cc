// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @detailed responsibilities:
/// @author Oliver Lange

// Unit Headers
#include <protocols/environment/claims/EnvClaim.hh>

// Package Headers
#include <protocols/environment/ClaimingMover.hh>
#include <core/environment/SequenceAnnotation.hh>

#include <protocols/environment/claims/ClaimStrength.hh>

// Project Headers

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

//// C++ headers
#include <iostream>
#include <utility/vector1.hh>

// option key includes


static basic::Tracer tr("protocols.environment.EnvClaim",basic::t_info);

namespace protocols {
namespace environment {
namespace claims {

EnvClaim::EnvClaim( ClaimingMoverOP owner ):
  ReferenceCount(),
  claim_source_( owner )
{}

/// @details Auto-generated virtual destructor
EnvClaim::~EnvClaim() {}

void EnvClaim::show( std::ostream& os ) const {
    os << "owned by, " << owner()->type() << ";";
}

ClaimingMoverOP EnvClaim::owner() const {
  return claim_source_;
}

extern std::ostream& operator<<( std::ostream& os, EnvClaim const& dof ) {
  dof.show( os );
  return os;
}

extern std::ostream& operator<<( std::ostream& os, EnvClaims const& dofs ) {
  for ( EnvClaims::const_iterator it = dofs.begin(); it != dofs.end(); ++it ) {
    if ( *it ) {
      os << **it << "\n";
    } else {
      os << "No-Claim\n";
    }
  }
  return os;
}

} //claims
} //environment
} //protocols
