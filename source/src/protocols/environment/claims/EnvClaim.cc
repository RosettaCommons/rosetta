// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/claims/EnvClaim.cc
/// @brief  Abstract class for implementing an EnvClaim (e.g. TorsionClaim).
/// @author Justin R. Porter

// Unit Headers
#include <protocols/environment/claims/EnvClaim.hh>

// Package Headers
#include <protocols/environment/ClaimingMover.hh>
#include <core/environment/SequenceAnnotation.hh>

#include <protocols/environment/claims/ClaimStrength.hh>

#include <protocols/environment/claims/CutBiasClaim.hh>
#include <protocols/environment/claims/JumpClaim.hh>
#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/XYZClaim.hh>

// Project Headers

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

#include <boost/algorithm/string.hpp>

//// C++ headers
#include <iostream>

// option key includes


static basic::Tracer tr("protocols.environment.EnvClaim",basic::t_info);

namespace protocols {
namespace environment {
namespace claims {

EnvClaimOP EnvClaim::make_claim( std::string const& name,
                                 ClaimingMoverOP owner,
                                 utility::tag::TagCOP tag,
                                 basic::datacache::DataMap& datamap ) {
  if      ( name == "CutBiasClaim" ) return new CutBiasClaim( owner, tag );
  else if ( name == "JumpClaim" )    return new JumpClaim( owner, tag );
  else if ( name == "TorsionClaim" ) return new TorsionClaim( owner, tag, datamap );
  else if ( name == "XYZClaim" )     return new XYZClaim( owner, tag, datamap );
  else throw utility::excn::EXCN_RosettaScriptsOption( "'" + name + "' is not a known EnvClaim type." );
}

bool EnvClaim::is_claim( std::string const& name ) {
  if      ( name == "CutBiasClaim" ) return true;
  else if ( name == "JumpClaim" )    return true;
  else if ( name == "TorsionClaim" ) return true;
  else if ( name == "XYZClaim" )     return true;
  else return false;
}

EnvClaim::EnvClaim( ClaimingMoverOP owner ):
  ReferenceCount(),
  claim_source_( owner )
{}

/// @details Auto-generated virtual destructor
EnvClaim::~EnvClaim() {}

void EnvClaim::show( std::ostream& os ) const {
    os << "owned by, " << owner()->type() << ";";
}

DOFElement EnvClaim::wrap_dof_id( core::id::DOF_ID const& id ) const {
  DOFElement e;
  e.id = id;

  return e;
}

ClaimingMoverOP EnvClaim::owner() const {
  return claim_source_;
}

ControlStrength EnvClaim::parse_ctrl_str( utility::tag::TagCOP tag ) const {
  return parse_ctrl_str( tag->getOption< std::string>( "control_strength" ) );
}

ControlStrength EnvClaim::parse_ctrl_str( std::string const& str ) const {
  std::string lower = str;
  boost::algorithm::to_lower( lower );

  if( lower == "does_not_control" ){
    return DOES_NOT_CONTROL;
  } else if( lower == "can_control" ){
    return CAN_CONTROL;
  } else if( lower == "must_control" ){
    return MUST_CONTROL;
  } else if( lower == "exclusive" ){
    return EXCLUSIVE;
  } else {
    throw utility::excn::EXCN_BadInput( "The initialization strength '" + str +
                                        "' is not recognized." );
  }
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
