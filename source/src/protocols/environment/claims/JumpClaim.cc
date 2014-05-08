// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file JumpClaim
/// @brief Claims access to a torsional angle.
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/claims/JumpClaim.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>
#include <core/environment/DofPassport.hh>

#include <protocols/environment/claims/EnvClaim.hh>

#include <protocols/environment/ClaimingMover.hh>
#include <protocols/environment/ProtectedConformation.hh>

// Project Headers
#include <core/id/TorsionID.hh>

#include <core/pose/util.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

// C++ headers

// option key includes

static basic::Tracer tr("protocols.environment.JumpClaim",basic::t_info);

namespace protocols {
namespace environment {
namespace claims {

using core::environment::LocalPosition;

JumpClaim::JumpClaim( ClaimingMoverOP owner,
                      std::string const& jump_label,
                      LocalPosition const& jpos1,
                      LocalPosition const& jpos2,
                      LocalPosition const& cutp ):
  EnvClaim( owner ),
  label_( jump_label ),
  pos1_( jpos1 ),
  pos2_( jpos2 ),
  cut_( cutp ),
  atom1_( "CA" ),
  atom2_( "CA" ),
  physical_cut_( false ),
  c_str_( EXCLUSIVE ),
  i_str_( DOES_NOT_INITIALIZE )
{}

void JumpClaim::yield_elements( core::environment::FoldTreeSketch const&, JumpElements& elements ) const {
  JumpElement e;

  e.p1 = pos1();
  e.p2 = pos2();
  e.label = label();
  e.atom1 = atom1();
  e.atom2 = atom2();
  e.has_physical_cut = physical();

  elements.push_back( e );
}

void JumpClaim::yield_elements( ProtectedConformationCOP const& conf, DOFElements& elements ) const {
  for( core::Size rb_i = core::id::RB1; rb_i <= core::id::RB6; ++rb_i ){
    DOFElement e;

    e.id = core::id::DOF_ID( conf->jump_atom_id( (int) conf->annotations()->resolve_jump( label() ) ),
                             core::id::DOF_Type( rb_i ) );

    e.i_str = i_str_;
    e.c_str = c_str_;
    e.owner.operator=( owner() );

    elements.push_back( e );
  }
}

void JumpClaim::yield_elements( core::environment::FoldTreeSketch const&, CutElements& elements ) const {
  if( cut_position() != LocalPosition( "", 0 ) ){
    CutElement e;
    e.p = cut_position();
    e.physical = physical();
    elements.push_back( e );
  } else if( physical() == true ){
    utility_exit_with_message("[FATAL] Physical (i.e. non-chainbreak-scored) jump claims are only avaliable explicitly specified cut positions.");
  }
}

void JumpClaim::set_atoms( std::string const& a1, std::string const& a2 ) {
  atom1_ = a1;
  atom2_ = a2;
}

void JumpClaim::ctrl_strength( ControlStrength const& str ){
  c_str_ = str;
}

void JumpClaim::init_strength( InitializationStrength const& str ){
  i_str_ = str;
}

std::string const& JumpClaim::label() const {
  return label_;
}

LocalPosition const& JumpClaim::pos1() const {
  return pos1_;
}

LocalPosition const& JumpClaim::pos2() const {
  return pos2_;
}

LocalPosition const& JumpClaim::cut_position() const {
  return cut_;
}

std::string const& JumpClaim::atom1() const {
  return atom1_;
}

std::string const& JumpClaim::atom2() const {
  return atom2_;
}

EnvClaimOP JumpClaim::clone() const {
  return new JumpClaim( *this );
}

std::string JumpClaim::str_type() const{
  return "Jump";
}

void JumpClaim::show( std::ostream& os ) const {
  os << str_type() << "(" << pos1() << "," << atom1() << "->" << pos2()
     << "," << atom2() << ") owned by a " << owner()->get_name();
}

} //claims
} //environment
} //protocols
