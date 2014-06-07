// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ClaimStrength.cc
/// @brief ClaimStrength tracks the prioirty of a particular claim with respect to its initialization
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/claims/ClaimStrength.hh>

// Package Headers

// Utility Headers
#include <core/types.hh>

// C++ Headers


namespace protocols {
namespace environment {
namespace claims {

ClaimStrength::ClaimStrength( PrioSubtype subtype, Size subprio ):
  Parent(),
  subtype_( subtype ),
  subprio_( subprio )
{}


ClaimStrength::ClaimStrength( ClaimStrength const& src ) :
  Parent(),
  subtype_( src.subtype_ ),
  subprio_( src.subprio_ )
{}


bool ClaimStrength::operator< ( ClaimStrength const& other ) const {
  if( other.subtype_ < subtype_ ){
    return true;
  } else if ( other.subtype_ > subtype_ ){
    return false;
  } else {
    if( other.subprio_ < subprio_ ){
      return true;
    } else {
      return false;
    }
  }
}


bool ClaimStrength::operator== ( ClaimStrength const& other ) const {
  if( other.subtype_ == subtype_ &&
      other.subprio_ == subprio_ ){
    return true;
  } else {
    return false;
  }
}

bool ClaimStrength::operator!= ( ClaimStrength const& other ) const{
  return !operator==( other );
}

bool ClaimStrength::operator> (ClaimStrength const& other ) const {
  return other.operator<( *this );
}

bool ClaimStrength::operator>= (ClaimStrength const& other) const {
  return !operator<( other );
}

bool ClaimStrength::operator<= ( ClaimStrength const& other ) const{
  return !operator>( other );
}

ClaimStrength::PrioSubtype ClaimStrength::subtype() const {
  return subtype_;
}

core::Size ClaimStrength::subprio() const {
  return subprio_;
}

extern std::ostream& operator<<( std::ostream& os, ClaimStrength const& ir) {
  os << "(" << ir.subtype() << "," << ir.subprio() << ")";
  return os;
}

extern std::ostream& operator<<( std::ostream& os, ClaimStrength::PrioSubtype const& ir) {
  if( ir == ClaimStrength::CAN_CONTROL ){
    os << "CAN_INIT";
  } else if ( ir == ClaimStrength::MUST_CONTROL ){
    os << "MUST_INIT";
  } else if ( ir == ClaimStrength::EXCLUSIVE ){
    os << "EXCLUSIVE";
  } else {
    assert( false );
  }
  return os;
}

} // protocols
} // environment
} // claims
