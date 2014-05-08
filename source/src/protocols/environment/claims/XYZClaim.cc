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
  c_str_( MUST_CONTROL),
  i_str_( DOES_NOT_INITIALIZE )
{}

XYZClaim::XYZClaim( ClaimingMoverOP owner,
                    LocalPosition const& local_pos):
  EnvClaim( owner ),
  c_str_( MUST_CONTROL ),
  i_str_( DOES_NOT_INITIALIZE )
{
  local_positions_ = LocalPositions();
  local_positions_.push_back( new LocalPosition( local_pos ) );
}

XYZClaim::XYZClaim( ClaimingMoverOP owner,
                            std::string const& label,
                            std::pair< core::Size, core::Size > const& range ):
  EnvClaim( owner ),
  c_str_( MUST_CONTROL ),
  i_str_( DOES_NOT_INITIALIZE )
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

void XYZClaim::build_bond_length_elements( core::Size seqpos,
                                           ProtectedConformationCOP const& conf,
                                           DOFElements& elements ) const {
  using namespace core::id;
  using namespace core::kinematics::tree;

  for( core::Size i = 1; i <= conf->residue( seqpos ).atoms().size(); ++i ){
    AtomID const a_id( i, seqpos );

    // This atom's not a BondedAtom, then this atom isn't
    // built by a bond (rather a jump) and doesn't have a well-defined bond length.
		core::kinematics::tree::Atom const& atom = conf->atom_tree().atom( a_id );

    if( !atom.is_jump() ){
      bool const is_external = ( atom.id().rsd() != atom.parent()->id().rsd() );

      if( claim_external() || !is_external ){
        DOF_ID id( a_id, core::id::D );
        elements.push_back( wrap_dof_id( id ) );
      }

      // If we're doing external bond lengths, one is length is "owned" by
      // the first atom in the next residue as well.
      if( claim_external() ){
        for( Size j = 0; j < atom.n_children(); ++j ){
          AtomCOP child = atom.child( j );
          if( child->id().rsd() != atom.id().rsd() ){
            elements.push_back( wrap_dof_id( DOF_ID( child->id(), core::id::D ) ) );
          }
        }
      }
    }
  }
}

void XYZClaim::build_bond_angle_elements( core::Size seqpos,
                                          ProtectedConformationCOP const& conf,
                                          DOFElements& elements ) const {
  using namespace core::id;
  using namespace core::kinematics::tree;

  for( Size i = 1; i <= conf->residue( seqpos ).atoms().size(); ++i ){
    AtomID const a_id( i, seqpos );
		core::kinematics::tree::Atom const& atom = conf->atom_tree().atom( a_id );

    if( !atom.is_jump() ){
      AtomCOP parent = atom.parent();
      if( parent && !parent->is_jump() && parent->parent() ){
        bool const angle_is_external = ( a_id.rsd()         != parent->id().rsd()           ) ||
                                       ( parent->id().rsd() != parent->parent()->id().rsd() );

        if( claim_external() || !angle_is_external ){
          // If we're doing external angles, always add. If we're not, add only if
          // the defined angle is fully internal.

          elements.push_back( wrap_dof_id( DOF_ID( a_id, core::id::THETA ) ) );
        }
      }
    }

    if( claim_external() ){
      // The position of the "last" atom in residue i helps determine the angles defined by
      // it's children and it's children's children. If we're doing external angles, do those too.
      for( Size j = 0; j < atom.n_children(); ++j ){
        AtomCOP child = atom.child( j );

        if( child->atom_id().rsd() != a_id.rsd() &&
            !child->is_jump() ){

          elements.push_back( wrap_dof_id(  DOF_ID( child->id(), core::id::THETA ) ) );

          for( Size k = 0; k < child->n_children(); ++k ){
            AtomCOP subchild = child->child(k);
            if( !subchild->is_jump() ){
              elements.push_back( wrap_dof_id( DOF_ID( subchild->id(), core::id::THETA ) ) );
            }
          }
        }
      }
    }
  }
}

void XYZClaim::build_bond_torsion_elements( core::Size seqpos,
                                            ProtectedConformationCOP const& conf,
                                            DOFElements& elements ) const {
  using namespace core::id;
  using namespace core::kinematics::tree;

  for( Size i = 1; i <= conf->residue( seqpos ).atoms().size(); ++i ){
    AtomID const a_id( i, seqpos );
		core::kinematics::tree::Atom const& a1 = conf->atom_tree().atom( a_id );
    std::cout << a_id << ":" << conf->residue( seqpos ).atom_name( i ) << " p@" << a1.parent() << std::endl;

    if( a1.parent() && !a1.parent()->is_jump() ){
      AtomCOP a2 = a1.parent();
      std::cout << "  " << a2->id() << ":" << conf->residue( seqpos ).atom_name( a2->id().atomno() ) << std::endl;

      if( a2->parent() && !a2->parent()->is_jump() ){
        AtomCOP a3 = a2->parent();
        std::cout << "  " << a3->id() << ":" << conf->residue( seqpos ).atom_name( a3->id().atomno() ) << std::endl;

        if( a3->parent() && !a3->parent()->is_jump() ){

          bool const angle_is_external = ( a1.id().rsd() != a2->id().rsd() ||
                                           a1.id().rsd() != a3->id().rsd() );
          if( claim_external() || !angle_is_external ){
            elements.push_back( wrap_dof_id( DOF_ID( a_id, core::id::PHI ) ) );
          }
        }
      }
    }

    if( claim_external() ){

      for( Size j = 0; j < a1.n_children(); ++j ){
        AtomCOP a2 = a1.child( j );
        if( a2->is_jump() ) continue;

        for( Size k = 0; k < a2->n_children(); ++k ){
          AtomCOP a3 = a2->child( k );
          if( a3->is_jump() ) continue;

          for( Size l = 0; l < a3->n_children(); ++l ){
            AtomCOP a4 = a3->child( l );
            if( a4->is_jump() ) continue;

            elements.push_back( wrap_dof_id( DOF_ID( a4->id(), core::id::PHI ) ) );
          }
        }
      }
    }
  }
}

void XYZClaim::yield_elements( ProtectedConformationCOP const& conf, DOFElements& elements ) const {

  for( LocalPositions::const_iterator pos = positions().begin();
      pos != positions().end(); ++pos ){

    Size seqpos = conf->annotations()->resolve_seq( **pos );
    for( Size i = 1; i <= conf->residue( seqpos ).atoms().size(); ++i ){
      elements.push_back( wrap_dof_id( core::id::DOF_ID( core::id::AtomID( i, seqpos ), core::id::D ) ) );
      elements.push_back( wrap_dof_id( core::id::DOF_ID( core::id::AtomID( i, seqpos ), core::id::THETA ) ) );
      elements.push_back( wrap_dof_id( core::id::DOF_ID( core::id::AtomID( i, seqpos ), core::id::PHI ) ) );
    }
  }
}

LocalPositions const& XYZClaim::positions() const {
  return local_positions_;
}

ControlStrength const& XYZClaim::ctrl_strength() const {
  return c_str_;
}

void XYZClaim::ctrl_strength( ControlStrength const& c_str ){
  if( c_str > EXCLUSIVE || c_str < DOES_NOT_CONTROL ){
    throw utility::excn::EXCN_RangeError( "ControlStrengths are limited to values between DOES_NOT_CONTROL and EXCLUSIVE" );
  } else {
    c_str_ = c_str;
  }
}

InitializationStrength const& XYZClaim::init_strength() const {
  return i_str_;
}

void XYZClaim::init_strength( InitializationStrength const& i_str ){
  if( i_str > MUST_INITIALIZE || i_str < DOES_NOT_INITIALIZE ){
    throw utility::excn::EXCN_RangeError( "InitializationStrengths are limited to values between DOES_NOT_INITIALIZE and MUST_INITIALIZE" );
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
