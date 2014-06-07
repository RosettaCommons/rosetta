// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TorsionClaim
/// @brief Claims access to a torsional angle.
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/claims/TorsionClaim.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>

#include <protocols/environment/claims/EnvClaim.hh>
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


static basic::Tracer tr("protocols.environment.TorsionClaim",basic::t_info);

namespace protocols {
namespace environment {
namespace claims {

using core::environment::LocalPosition;
using core::environment::LocalPositions;

TorsionClaim::TorsionClaim( ClaimingMoverOP owner,
                            LocalPosition const& local_pos):
  EnvClaim( owner ),
  c_str_( MUST_CONTROL ),
  i_str_( DOES_NOT_CONTROL ),
  claim_sidechain_( false ),
  claim_backbone_( true )
{
  local_positions_ = LocalPositions();
  local_positions_.push_back( new LocalPosition( local_pos ) );
}

TorsionClaim::TorsionClaim( ClaimingMoverOP owner,
                            std::string const& label,
                            std::pair< core::Size, core::Size > const& range ):
  EnvClaim( owner ),
  c_str_( MUST_CONTROL ),
  i_str_( DOES_NOT_CONTROL ),
  claim_sidechain_( false ),
  claim_backbone_( true )
{
  local_positions_ = LocalPositions();
  for( Size i = range.first; i <= range.second; ++i){
    local_positions_.push_back( new LocalPosition( label, i ) );
  }
}

TorsionClaim::TorsionClaim( ClaimingMoverOP owner,
                           LocalPositions const& positions ):
  EnvClaim( owner ),
  local_positions_( positions ),
  c_str_( MUST_CONTROL ),
  i_str_( DOES_NOT_CONTROL ),
  claim_sidechain_( false ),
  claim_backbone_( true )
 {}


DOFElement TorsionClaim::wrap_dof_id( core::id::DOF_ID const& id ) const {
  DOFElement e = Parent::wrap_dof_id( id );

  e.c_str = ctrl_strength();
  e.i_str = init_strength();

  return e;
}

void TorsionClaim::yield_elements( ProtectedConformationCOP const& conf, DOFElements& elements ) const {

  using namespace core::id;

  for( LocalPositions::const_iterator pos = positions().begin();
       pos != positions().end(); ++pos ){

    Size seqpos = conf->annotations()->resolve_seq( **pos );

    if( conf->residue( seqpos ).type().is_protein() ) {
      if( claim_backbone() ){
        TorsionID phi  ( seqpos, BB, phi_torsion );
        TorsionID psi  ( seqpos, BB, psi_torsion );
        TorsionID omega( seqpos, BB, omega_torsion );

        DOF_ID phi_dof = conf->dof_id_from_torsion_id( phi );
        DOF_ID psi_dof = conf->dof_id_from_torsion_id( psi );
        DOF_ID omega_dof = conf->dof_id_from_torsion_id( omega );

        elements.push_back( wrap_dof_id( phi_dof ) );
        elements.push_back( wrap_dof_id( psi_dof ) );
        elements.push_back( wrap_dof_id( omega_dof ) );
      }
      if( claim_sidechain() ){
        for( Size i = 1; i <= conf->residue( seqpos ).nchi(); ++i ){
          TorsionID t_id( seqpos, CHI, i );
          elements.push_back( wrap_dof_id( conf->dof_id_from_torsion_id( t_id ) ) );
        }
      }
    } else {
      // Hey there! If you got here because you're using sugars or RNA or something,
      // all you have to do is put an else if for your case, and correctly
      // generate the DOF_IDs for your case and then use wrap_dof_id to build
      // a DOFElement out of it. JRP
      tr.Debug << "TorsionClaim owned by " << owner()->get_name()
               << " is ignoring residue " << seqpos
               << " with name3 " << conf->residue( seqpos ).name3()
               << " because it is not protein." << std::endl;
    }
  }
}


LocalPositions const& TorsionClaim::positions() const {
  return local_positions_;
}
void TorsionClaim::strength( ControlStrength const& c_str, ControlStrength const& i_str ){
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

ControlStrength const& TorsionClaim::ctrl_strength() const {
  return c_str_;
}


ControlStrength const& TorsionClaim::init_strength() const {
  return i_str_;
}

EnvClaimOP TorsionClaim::clone() const {
  return new TorsionClaim( *this );
}

std::string TorsionClaim::str_type() const{
  return "Torsion";
}

void TorsionClaim::show( std::ostream& os ) const {
  os << str_type() << " owned by a " << owner()->get_name() << " with " << positions().size() << "positions.";
}

} //claims
} //environment
} //protocols
