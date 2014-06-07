// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ClaimStrength.hh
/// @brief Header file for ClaimStrength, which tracks the prioirty of a particular claim with respect to its initialization
/// @author Justin Porter


#ifndef INCLUDED_protocols_environment_claims_ClaimStrength_hh
#define INCLUDED_protocols_environment_claims_ClaimStrength_hh


// Unit Headers
#include <protocols/environment/claims/ClaimStrength.fwd.hh>

// Package Headers

// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <string>
#include <ostream>

namespace protocols {
namespace environment {
namespace claims {

class ClaimStrength : public utility::pointer::ReferenceCount {
  typedef utility::pointer::ReferenceCount Parent;

public:
  // PrioSubtypes indicate the level of control over the dof desired by the mover.
  // DOES_NOT_CONTROL is the lowest level, used for protocols that exert no control
  //   over the claimed dof (perhaps, for example, because they only require that it exists).
  // CAN_CONTROL is also a low, uncommon level indicating the mover would like to control
  //   (i.e. samples) the dof, but needn't. Fragment Insertion is often reliant on claims
  //   of this type, because other movers can take priority (i.e. if a region is fixed).
  // MUST_CONTROL asserts a need to sample the degree of freedom (e.g. jumps in docking).
  //   This level allows the mover, however, to share the dof with other movers.
  // EXCLUSIVE is as MUST_CONTROL, except that no other mover can be granted access to the dof.
  enum PrioSubtype  {
    DOES_NOT_CONTROL = 0,
    CAN_CONTROL = 1,
    MUST_CONTROL,
    EXCLUSIVE
  };

  ClaimStrength( PrioSubtype type, Size subprio = 0 );

  ClaimStrength( ClaimStrength const& );

  bool operator<( ClaimStrength const& ) const;

  bool operator==( ClaimStrength const& ) const;

  bool operator!=( ClaimStrength const& ) const;

  bool operator>( ClaimStrength const& ) const;

  bool operator>= (ClaimStrength const& ) const;

  bool operator<= ( ClaimStrength const& ) const;

  PrioSubtype subtype() const;

  Size subprio() const;

private:

  PrioSubtype subtype_;
  Size subprio_;

};

extern std::ostream& operator<<( std::ostream& os, ClaimStrength const& ir );
extern std::ostream& operator<<( std::ostream& os, ClaimStrength::PrioSubtype const& ir );


}
}
}

#endif
