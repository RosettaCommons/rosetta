// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/DofPassport.cc
/// @author Justin R. Porter

// Unit Headers
#include <core/environment/DofPassport.hh>

// Package headers
#include <core/environment/EnvCore.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/Conformation.hh>

#include <utility/excn/Exceptions.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr("protocols.environment.DofPassport", basic::t_info);

namespace core {
namespace environment {

using core::kinematics::MoveMapOP;

DofPassport::DofPassport( std::string const& mover,
                         core::Size env ) :
  mover_( mover ),
  env_id_( env )
{}

DofPassport::~DofPassport() {}

std::set< core::id::DOF_ID >::const_iterator DofPassport::begin() const {
  if( accessible_dofs_.begin()->rsd() == 0 ){
    return ++accessible_dofs_.begin();
  } else {
    return accessible_dofs_.begin();
  }

}


//Movemap Config
MoveMapOP DofPassport::render_movemap() const {
  MoveMapOP mm = new core::kinematics::MoveMap();

  for( std::set<core::id::DOF_ID>::iterator id_it = accessible_dofs_.begin();
       id_it != accessible_dofs_.end(); ++id_it ){
    mm->set( *id_it, true);
  }

  for( core::Size seqpos = 1; seqpos <= conf_->size(); ++seqpos ){
    if( conf_->residue( seqpos ).is_protein() ){
      bool seqpos_access = true;
      for( Size torsion_type = id::phi_torsion; torsion_type <= id::omega_torsion; ++torsion_type ){
        id::TorsionID t_id = id::TorsionID( seqpos, id::BB, torsion_type );
        id::DOF_ID d_id = conf_->dof_id_from_torsion_id( t_id );

        seqpos_access &= d_id.valid() && accessible_dofs_.find( d_id ) != accessible_dofs_.end();
      }
      mm->set_bb( seqpos, seqpos_access );

      seqpos_access = true;
      for( core::Size chi_i = 1; chi_i <= conf_->residue( seqpos ).nchi(); ++chi_i ){
        id::TorsionID t_id = id::TorsionID( seqpos, id::CHI, chi_i );
        id::DOF_ID d_id = conf_->dof_id_from_torsion_id( t_id );

        seqpos_access &= d_id.valid() && accessible_dofs_.find( d_id ) != accessible_dofs_.end();
      }
      mm->set_chi( seqpos, seqpos_access );
    } else { // nonprotein (e.g. VRT) can't be configured into a movemap?
      // TODO: add cases for sugars, RNA, etc.
      mm->set_bb( seqpos, false );
      mm->set_chi( seqpos, false );
      if( !conf_->residue( seqpos ).is_virtual_residue() &&
          !conf_->residue( seqpos ).is_virtual( 1 ) ){
        tr.Warning << "Residue " << seqpos << " named " << conf_->residue( seqpos ).name3()
                   << " is being ignored by DofPassport::render (" << __FILE__ << ":"
                   << __LINE__ << ") because it doesn't how to encode information about its "
                   << "backbone and sidechain angles in the movemap." << std::endl;
      }
    }
  }

  for( int jump_i = 1; jump_i <= (int) conf_->fold_tree().num_jump(); ++jump_i ){
    bool allow = has_jump_access( jump_i );
    mm->set_jump( jump_i, allow );
    mm->set_jump( conf_->fold_tree().upstream_jump_residue( jump_i ),
                  conf_->fold_tree().downstream_jump_residue( jump_i ),
                  allow );
  }

  return mm;
}

bool DofPassport::has_jump_access( int jump_num ) const {
  bool allow = true;
  for( core::Size rb_i = id::RB1; rb_i <= id::RB6; ++rb_i ){
    id::AtomID a_id = conf_->jump_atom_id( jump_num );
    allow &= accessible_dofs_.find( id::DOF_ID( a_id, id::DOF_Type(rb_i) ) ) != accessible_dofs_.end();
  }
  return allow;
}

utility::vector1< int > DofPassport::active_jumps() const {
  utility::vector1< int > active_jumps;
  for( int i = 1; i <= (int) conf_->fold_tree().num_jump(); ++i ){
    if( has_jump_access( i ) ){
      active_jumps.push_back( i );
    }
  }
  return active_jumps;
}

void DofPassport::add_dof_access( id::DOF_ID const& dof_id ){
  accessible_dofs_.insert( dof_id );
}

void DofPassport::revoke_all_access(){
  accessible_dofs_.clear();
}

//Security Checks
bool DofPassport::access_check( EnvCore const& env, bool type_specific_check ) const {
  if( env.id() == env_id_ ){
    if( type_specific_check ){
      if( tr.Trace.visible() ) tr.Trace << "allowed" << std::endl;
      return true;
    } else {
      if( tr.Trace.visible() ) tr.Trace << "denied (right not given)" << std::endl;
      return false;
    }
  } else {
    if( tr.Trace.visible() ) tr.Trace << "denied (environmental conflict)" << std::endl;
    return false;
  }
}
bool DofPassport::dof_access( id::DOF_ID const& id ) const {
  return ( accessible_dofs_.find( id ) != accessible_dofs_.end() );
}

bool DofPassport::dof_access( EnvCore const& env, id::DOF_ID const& id ) const {
  if( tr.Trace.visible() ) tr.Trace << "Verifying access to DOF_ID " << id << ": ";
  return access_check( env, dof_access( id ) );
}

void DofPassport::reference_conformation( conformation::ConformationCOP conf ) {
  conf_ = conf;
}

conformation::ConformationCOP DofPassport::reference_conformation() const {
  return conf_;
}

std::string const& DofPassport::mover() const {
  return mover_;
}

void DofPassport::show( std::ostream& str ) const {
  str << "DofPassport: " << mover_ << ", " << env_id_;
}

std::ostream& operator<<( std::ostream& str, DofPassport const& dof_passport ) {
  dof_passport.show( str );
  return str;
}

} // environment
} // core
