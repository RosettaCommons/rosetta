// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file XYZClaim
/// @brief Claims access to a torsional angle.
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/claims/XYZClaim.hh>

// Package Headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>

#include <core/environment/LocalPosition.hh>
#include <core/environment/SequenceAnnotation.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/ClaimStrength.hh>

#include <protocols/environment/ClaimingMover.hh>
#include <protocols/environment/ProtectedConformation.hh>

// Project Headers
#include <core/id/TorsionID.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

// C++ headers

// option key includes


static basic::Tracer tr("protocols.environment.XYZClaim",basic::t_info);

namespace protocols {
namespace environment {
namespace claims {

using core::environment::LocalPosition;
using core::environment::LocalPositions;

XYZClaim::XYZClaim( ClaimingMoverOP owner ):
  EnvClaim( owner ),
  c_str_( DOES_NOT_CONTROL),
  i_str_( DOES_NOT_CONTROL )
{}

XYZClaim::XYZClaim( ClaimingMoverOP owner,
                    LocalPosition const& local_pos):
  EnvClaim( owner ),
  c_str_( DOES_NOT_CONTROL ),
  i_str_( DOES_NOT_CONTROL )
{
  local_positions_ = LocalPositions();
  local_positions_.push_back( new LocalPosition( local_pos ) );
}

XYZClaim::XYZClaim( ClaimingMoverOP owner,
                            std::string const& label,
                            std::pair< core::Size, core::Size > const& range ):
  EnvClaim( owner ),
  c_str_( DOES_NOT_CONTROL ),
  i_str_( DOES_NOT_CONTROL )
{
  local_positions_ = LocalPositions();
  for( Size i = range.first; i <= range.second; ++i){
    local_positions_.push_back( new LocalPosition( label, i ) );
  }
}

void XYZClaim::add_position( LocalPosition const& p ){
  local_positions_.push_back( new LocalPosition( p ) );
}


DOFElement XYZClaim::wrap_dof_id( core::id::DOF_ID const& id ) const {
  DOFElement e = Parent::wrap_dof_id( id );

  e.c_str = ctrl_strength();
  e.i_str = init_strength();

  return e;
}

void XYZClaim::yield_elements( ProtectedConformationCOP const& conf, DOFElements& elements ) const {

  if( ctrl_strength() == DOES_NOT_CONTROL &&
      init_strength() == DOES_NOT_CONTROL ){
    tr.Warning << "XYZClaim owned by " << owner()->get_name()
               << " has both initializaiton and sampling strength set to DOES_NOT_CONTROL."
               << "  Did you forget to set these values?" << std::endl;
  }

  for( LocalPositions::const_iterator pos = positions().begin();
      pos != positions().end(); ++pos ){

    Size seqpos = conf->annotations()->resolve_seq( **pos );
    for( Size i = 1; i <= conf->residue( seqpos ).atoms().size(); ++i ){
      core::id::AtomID atom_id( i, seqpos );
      if( conf->atom_tree().atom( atom_id ).is_jump() ){
        for( int j = core::id::RB1; j <= core::id::RB6; ++j ){
          elements.push_back( wrap_dof_id( core::id::DOF_ID( atom_id,
                                                             core::id::DOF_Type( j ) ) ) );
        }
      } else {
        for( int j = core::id::D; j <= core::id::PHI; ++j ){
          elements.push_back( wrap_dof_id( core::id::DOF_ID( atom_id,
                                                             core::id::DOF_Type( j ) ) ) );
        }
      }
    }
  }
}

LocalPositions const& XYZClaim::positions() const {
  return local_positions_;
}

ControlStrength const& XYZClaim::ctrl_strength() const {
  return c_str_;
}

ControlStrength const& XYZClaim::init_strength() const {
  return i_str_;
}

void XYZClaim::strength( ControlStrength const& c_str, ControlStrength const& i_str ){
  if( c_str > EXCLUSIVE || c_str < DOES_NOT_CONTROL ){
    throw utility::excn::EXCN_RangeError( "Sampling ControlStrengths are limited to values between DOES_NOT_CONTROL and EXCLUSIVE" );
  } else {
    c_str_ = c_str;
  }

  if( i_str > EXCLUSIVE || i_str < DOES_NOT_CONTROL ){
    throw utility::excn::EXCN_RangeError( "Initialization ControlStrengths are limited to values between DOES_NOT_INITIALIZE and EXCLUSIVE" );
  } else {
    i_str_ = i_str;
  }
}

EnvClaimOP XYZClaim::clone() const {
  return new XYZClaim( *this );
}

std::string XYZClaim::str_type() const{
  return "Torsion";
}

void XYZClaim::show( std::ostream& os ) const {
  os << str_type() << " owned by a " << owner()->get_name() << " with " << positions().size() << "positions.";
}

} //claims
} //environment
} //protocols
