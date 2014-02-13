// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/DofPassport.cc
/// @author Justin Porter

// Unit Headers
#include <core/environment/DofPassport.hh>

// Package headers
#include <core/environment/EnvCore.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>
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

//Movemap Config
MoveMapOP DofPassport::render_movemap() const{
  MoveMapOP mm = new core::kinematics::MoveMap();

  for( std::set<core::id::DOF_ID>::iterator id_it = accessible_dofs_.begin();
       id_it != accessible_dofs_.end(); ++id_it ){
    mm->set( *id_it, true);
    if( id_it->type() == core::id::RB1 ){
      mm->set( core::id::DOF_ID( id_it->atom_id(), core::id::RB2 ), true );
      mm->set( core::id::DOF_ID( id_it->atom_id(), core::id::RB3 ), true );
      mm->set( core::id::DOF_ID( id_it->atom_id(), core::id::RB4 ), true );
      mm->set( core::id::DOF_ID( id_it->atom_id(), core::id::RB5 ), true );
      mm->set( core::id::DOF_ID( id_it->atom_id(), core::id::RB6 ), true );
    }
  }

  for( std::set< id::TorsionID >::iterator id_it = accessible_torsions_.begin();
       id_it != accessible_torsions_.end(); ++id_it ){
    //assumption: if one torsion is accessible, all are accessible
    mm->set_bb( id_it->rsd(), true );
  }

  for( std::set<core::Size>::iterator id_it = accessible_jump_numbers_.begin();
       id_it != accessible_jump_numbers_.end(); ++id_it ){
    mm->set_jump( (int) *id_it, true );
  }

  for( std::set< id::JumpID >::const_iterator id_it = accessible_jump_ids_.begin();
      id_it != accessible_jump_ids_.end(); ++id_it ){
    mm->set_jump( *id_it, true );
  }

  return mm;
}

//Access right setup

void DofPassport::add_jump_access( id::AtomID const& id, Size const& nr, id::JumpID const& jid ){
  //Only RB1 is stored because I can't think of any reason RB2-6 don't follow RB1.
  accessible_dofs_.insert( id::DOF_ID( id, id::RB1 ) );
  accessible_jump_numbers_.insert( nr );
  accessible_jump_ids_.insert( jid );
}

void DofPassport::add_dof_access( id::DOF_ID const& dof_id, core::id::TorsionID const& torsion_id ){
  assert( dof_id.valid() && torsion_id.valid() );
  if( dof_id.type() >= id::RB1 && dof_id.type() <= id::RB6 ){
    throw utility::excn::EXCN_BadInput( "DofPassport being asked to add a jump through the add_dof_access function. Use add_jump_access instead." );
  } else {
    accessible_dofs_.insert( dof_id );
    accessible_torsions_.insert( torsion_id );
  }
}

void DofPassport::revoke_all_access(){
  accessible_dofs_.clear();
  accessible_jump_numbers_.clear();
  accessible_torsions_.clear();
  accessible_jump_ids_.clear();
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

bool DofPassport::jump_access( EnvCore const& env, id::AtomID const& atom_id ) const {
  if( tr.Trace.visible() ) tr.Trace << "Verifying access to jump from atom " << atom_id << ": ";
  id::DOF_ID id = id::DOF_ID( atom_id, id::RB1 );
  bool grant = ( accessible_dofs_.find( id ) != accessible_dofs_.end() );
  return access_check( env, grant );
}

bool DofPassport::jump_access( EnvCore const& env, Size const& nr ) const {
  if( tr.Trace.visible() ) tr.Trace << "Verifying access to jump number " << nr << ": ";
  bool grant = ( accessible_jump_numbers_.find( nr ) != accessible_jump_numbers_.end() );
  return access_check( env, grant );
}

bool DofPassport::dof_access( EnvCore const& env, id::DOF_ID const& id ) const {
  if( tr.Trace.visible() ) tr.Trace << "Verifying access to DOF_ID " << id << ": ";
  bool grant = (accessible_dofs_.find( id ) != accessible_dofs_.end() );
  return access_check( env, grant );
}

bool DofPassport::jump_access( EnvCore const& env, id::JumpID const& id ) const{
  if( tr.Trace.visible() ) tr.Trace << "Verifying access to JumpID " << id << ": ";
  bool grant = (accessible_jump_ids_.find( id ) != accessible_jump_ids_.end() );
  return access_check( env, grant );
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
