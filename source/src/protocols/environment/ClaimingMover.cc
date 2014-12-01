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
#include <protocols/environment/DofUnlock.hh>
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

static thread_local basic::Tracer tr( "protocols.environment.ClaimingMover", basic::t_info );

namespace protocols {
namespace environment {

using core::environment::DofPassport;
using core::environment::DofPassportCOP;

core::Real const TOLERANCE = 1e-6;

ClaimingMover::ClaimingMover():
  Mover(),
  passports_()
{}

ClaimingMover::ClaimingMover( ClaimingMover const& other ):
  Mover( other ),
  passports_( other.passports_ )
{
}

ClaimingMover::~ClaimingMover() {}

bool ang_delta( core::Real const& a, core::Real const& b ){
  return std::abs( std::cos( a ) - std::cos( b ) ) > TOLERANCE;
}

void ClaimingMover::sandboxed_copy( core::pose::Pose const& sandbox_pose,
                                    core::pose::Pose& true_pose ) const {
  // Copy result into protected conformation in in_pose
  environment::DofUnlock unlock( true_pose.conformation(), passport() );

  assert( sandbox_pose.conformation().fold_tree() == true_pose.conformation().fold_tree() );

  for ( Size i = 1; i <= true_pose.total_residue(); ++i ) {
    try {
      if( true_pose.residue( i ).is_protein() ){
        if( ang_delta( sandbox_pose.omega( i ), true_pose.omega( i ) ) ){
          true_pose.set_omega( i, sandbox_pose.omega( i ) );
        }
        if( ang_delta( sandbox_pose.phi( i ), true_pose.phi( i ) ) ){
          true_pose.set_phi( i, sandbox_pose.phi( i ) );
        }
        if( ang_delta( sandbox_pose.psi( i ), true_pose.psi( i ) ) ) {
          true_pose.set_psi( i, sandbox_pose.psi( i ) );
        }
        for( Size j = 1; j <= true_pose.conformation().residue( i ).nchi(); ++j ){
          if( ang_delta( true_pose.chi( (int) j, i ), sandbox_pose.chi( (int) j, i ) ) ){
            true_pose.set_chi( (int) j, i, sandbox_pose.chi( (int) j, i ) );
          }
        }
      }
    } catch ( environment::EXCN_Env_Security_Exception const& e ){
      tr.Error << "[ERROR] Unauthorized changes occurred during loop closure by mover '" << this->get_name()
      << "': (attempt to write to resid " << i << ")." << std::endl;
      throw e;
    }
  }

  for ( Size i = 1; i <= true_pose.num_jump(); ++i ){
    core::Real const delta_trans = sandbox_pose.jump( (int) i ).get_translation().distance( true_pose.jump( (int) i).get_translation() );
    numeric::xyzMatrix< double > const delta_rot = sandbox_pose.jump( (int) i ).get_rotation()*true_pose.jump( (int) i ).get_rotation().inverse();
    if( delta_trans < TOLERANCE && ( delta_rot.trace() - 3 ) < TOLERANCE ){
      try {
        true_pose.set_jump( (int) i, sandbox_pose.jump( (int) i ) );
      } catch ( environment::EXCN_Env_Security_Exception const& e ){
        tr.Error << "[ERROR] Unauthorized changes occurred during loop closure by mover '" << this->get_name()
        << "': (attempt to write to jump id " << i << ")." << std::endl;
        throw e;
      }
    }
  }
}


//CLAIMING METHODS:

void ClaimingMover::initialize( Pose& ) {
  throw utility::excn::EXCN_Msg_Exception( "ClaimingMover "+this->get_name()+
                                           " claimed a dof to initialize, but did not override ClaimingMover::initialize" );
}


//PASSPORT MANAGEMENT METHODS:
core::environment::DofPassportCOP ClaimingMover::passport() const {
  if( passports_.empty() ){
    throw utility::excn::EXCN_NullPointer( "ClaimingMover "+this->get_name()+
                                           " tried to access its passports, of which it has none.");
  }
  return passports_.top().second;
}

EnvironmentCAP ClaimingMover::active_environment() const {
  if( passports_.empty()) {
    return EnvironmentCAP();
  } else {
    return passports_.top().first;
  }
}


void ClaimingMover::push_passport( EnvironmentCAP env_ap, DofPassportCOP pass ){
  // Bad-state checks:
  // 1) If there's a superenvironment that mover is registered with, there SHOULD be a superenv passport
  // 2) If not, there should be an empty passport stack.
  EnvironmentCOP env( env_ap );
  EnvironmentCOP superenv( env->superenv().lock() );
  if( superenv && superenv->is_registered( utility::pointer::static_pointer_cast< ClaimingMover >( get_self_ptr() ) ) ){
		EnvironmentCOP passport_env = passports_.top().first.lock();
    if( passport_env && passport_env->id() == superenv->id() ){
      throw EXCN_Env_Passport( "ClaimingMover being double-assigned a passport for an environment.",
                               get_name(), env_ap );
    }
  } else {
    if( !( passports_.empty() ) ){
      if( passports_.top().first.lock()->id() == env->id() ) {
        std::ostringstream ss;
        ss << "ClaimingMover " << this->get_name() << " already has a passport for " << env->name()
           << ". This probably means Environment::start got called multiple times." << std::endl;
        throw EXCN_Env_Passport( ss.str(), get_name(), env_ap );
      } else {
        throw EXCN_Env_Passport( "ClaimingMover lacks a superenvironment passport for a superenvironment with which it is registered.",
                               get_name(), env_ap );
      }
    }
  }

  passports_.push( std::make_pair( env_ap, pass ) );
}

void ClaimingMover::pop_passport( EnvironmentCAP env_ap ) {
  // passport_cancel should never be called on an empty passport stack
  if( passports_.empty() ) {
    throw EXCN_Env_Passport( "Environment trying to pop a passport on an empty pass stack.", get_name(), env_ap );
  }
  // in fact, we should only ever cancel our own passports. (Should this be an assert?)
  EnvironmentCOP env( env_ap );
  EnvironmentCOP passport_env = passports_.top().first.lock();
  if( passport_env->id() != env->id() ){
    throw EXCN_Env_Passport( "Environment trying to pop a passport it did not issue.", get_name(), env_ap );
  }
  passports_.pop();
  passport_updated();
}

bool ClaimingMover::state_check( std::string const& method_name, bool test ) const {
  if( !has_passport() || test ){
    return true;
  } else {
    std::ostringstream ss;
    ss << "Call to " << this->get_name() << "::" << method_name << " is illegal because it has already "
    << "issued its claims." << std::endl;
    throw utility::excn::EXCN_Msg_Exception( ss.str() );
  }
}


} // environment
} // protocols
