// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/ClaimingMover.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/ClaimingMover.hh>

// Package headers
#include <core/environment/DofPassport.hh>
#include <protocols/environment/ProtectedConformation.hh>
#include <protocols/environment/Environment.hh>

#include <protocols/environment/EnvExcn.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Conformation.fwd.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr("protocols.environment.ClaimingMover", basic::t_info);

namespace protocols {
namespace environment {

using core::environment::DofPassport;
using core::environment::DofPassportCOP;

  ClaimingMover::ClaimingMover():
  Mover(),
  passports_()
{}

ClaimingMover::ClaimingMover( ClaimingMover const& other ):
  Mover( other ),
  passports_( other.passports_ )
{}

ClaimingMover::~ClaimingMover() {}

//CLAIMING METHODS:

void ClaimingMover::initialize( Pose& ) {
  throw utility::excn::EXCN_Msg_Exception( "ClaimingMover "+this->name()+
                                           " claimed a dof to initialize, but did not override ClaimingMover::initialize" );
}


//PASSPORT MANAGEMENT METHODS:
core::environment::DofPassportCOP ClaimingMover::passport() const {
  if( passports_.empty() ){
    throw utility::excn::EXCN_NullPointer( "ClaimingMover "+this->name()+
                                           " tried to access its passports, of which it has none.");
  }
  return passports_.top().second;
}

void ClaimingMover::push_passport( EnvironmentCAP env, DofPassportCOP pass ){
  // Bad-state checks:
  // 1) If there's a superenvironment that mover is registered with, there SHOULD be a superenv passport
  // 2) If not, there should be an empty passport stack.
  if( env->superenv() && env->superenv()->is_registered( this ) ){
    if( passports_.top().first == env->superenv()->id() ){
      throw EXCN_Env_Passport( "ClaimingMover being double-assigned a passport for an environment.",
                               get_name(), env );
    }
  } else {
    if( !( passports_.empty() ) ){
      throw EXCN_Env_Passport( "ClaimingMover lacks a superenvironment passport for a superenvironment with which it is registerd.",
                               get_name(), env );
    }
  }

  passports_.push( std::make_pair( env->id(), pass ) );
  passport_updated();
}

void ClaimingMover::pop_passport( EnvironmentCAP env ) {
  // passport_cancel should never be called on an empty passport stack
  if( passports_.empty() ) {
    throw EXCN_Env_Passport( "Environment trying to pop a passport on an empty pass stack.", get_name(), env );
  }
  // in fact, we should only ever cancel our own passports. (Should this be an assert?)
  if( passports_.top().first != env->id() ){
    throw EXCN_Env_Passport( "Environment trying to pop a passport it did not issue.", get_name(), env );
  }

  passports_.pop();
  passport_updated();
}

} // environment
} // protocols
