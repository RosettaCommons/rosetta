// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/environment/ClientMover.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/ClientMover.hh>

// Package headers
#include <core/environment/DofPassport.hh>
#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/ProtectedConformation.hh>
#include <protocols/environment/Environment.hh>

#include <protocols/environment/EnvExcn.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Conformation.fwd.hh>

#include <utility/excn/Exceptions.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr( "protocols.environment.ClientMover", basic::t_info );

namespace protocols {
namespace environment {

using core::environment::DofPassport;
using core::environment::DofPassportCOP;

core::Real const TOLERANCE = 1e-6;

ClientMover::ClientMover():
	Mover(),
	passports_()
{}

ClientMover::ClientMover( ClientMover const& )= default;

ClientMover::~ClientMover() = default;

bool ang_delta( core::Real const& a, core::Real const& b ){
	return std::abs( std::cos( a ) - std::cos( b ) ) > TOLERANCE;
}

void ClientMover::sandboxed_copy( core::pose::Pose const& sandbox_pose,
	core::pose::Pose& true_pose ) const {
	// Copy result into protected conformation in in_pose
	environment::DofUnlock unlock( true_pose.conformation(), passport() );

	debug_assert( sandbox_pose.conformation().fold_tree() == true_pose.conformation().fold_tree() );

	for ( Size i = 1; i <= true_pose.size(); ++i ) {
		try {
			if ( true_pose.residue( i ).is_protein() ) {
				if ( ang_delta( sandbox_pose.omega( i ), true_pose.omega( i ) ) ) {
					true_pose.set_omega( i, sandbox_pose.omega( i ) );
				}
				if ( ang_delta( sandbox_pose.phi( i ), true_pose.phi( i ) ) ) {
					true_pose.set_phi( i, sandbox_pose.phi( i ) );
				}
				if ( ang_delta( sandbox_pose.psi( i ), true_pose.psi( i ) ) ) {
					true_pose.set_psi( i, sandbox_pose.psi( i ) );
				}
				for ( Size j = 1; j <= true_pose.conformation().residue( i ).nchi(); ++j ) {
					if ( ang_delta( true_pose.chi( (int) j, i ), sandbox_pose.chi( (int) j, i ) ) ) {
						true_pose.set_chi( (int) j, i, sandbox_pose.chi( (int) j, i ) );
					}
				}
			}
		} catch ( environment::EXCN_Env_Security_Exception const& e ){
			tr.Error << "Unauthorized changes occurred during loop closure by mover '" << this->get_name()
				<< "': (attempt to write to resid " << i << ")." << std::endl;
			throw;// e;
		}
	}

	for ( Size i = 1; i <= true_pose.num_jump(); ++i ) {
		core::Real const delta_trans = sandbox_pose.jump( (int) i ).get_translation().distance( true_pose.jump( (int) i).get_translation() );
		numeric::xyzMatrix< double > const delta_rot = sandbox_pose.jump( (int) i ).get_rotation()*true_pose.jump( (int) i ).get_rotation().inverse();
		if ( delta_trans > TOLERANCE ||
				( delta_rot.trace() - 3 ) > TOLERANCE ) {
			try {
				true_pose.set_jump( (int) i, sandbox_pose.jump( (int) i ) );
			} catch ( environment::EXCN_Env_Security_Exception& e ){
				std::ostringstream ss;
				ss << "A differing value was detected in jump id " << i << " by " << __FUNCTION__ << ". "
					<< "Jump Translation difference: " << delta_trans << "; Rotation Difference: " << ( delta_rot.trace() - 3 )
					<< ".";
				e.add_msg( ss.str() );
				throw;// e;
			}
		}
	}
}


//CLAIMING METHODS:

void ClientMover::initialize( Pose& ) {
	throw CREATE_EXCEPTION(utility::excn::Exception,  "ClientMover "+this->get_name()+
		" claimed a dof to initialize, but did not override ClientMover::initialize" );
}


//PASSPORT MANAGEMENT METHODS:
core::environment::DofPassportCOP ClientMover::passport() const {
	if ( passports_.empty() ) {
		throw CREATE_EXCEPTION(utility::excn::NullPointerError,  "ClientMover "+this->get_name()+
			" tried to access its passports, of which it has none.");
	}
	return passports_.top().second;
}

EnvironmentCAP ClientMover::active_environment() const {
	if ( passports_.empty() ) {
		return EnvironmentCAP();
	} else {
		return passports_.top().first;
	}
}


void ClientMover::push_passport( EnvironmentCAP env_ap, DofPassportCOP pass ){
	// Bad-state checks:
	// 1) If there's a superenvironment that mover is registered with, there SHOULD be a superenv passport
	// 2) If not, there should be an empty passport stack.
	EnvironmentCOP env( env_ap );
	EnvironmentCOP superenv( env->superenv().lock() );
	if ( superenv && superenv->is_registered( utility::pointer::static_pointer_cast< ClientMover >( get_self_ptr() ) ) ) {
		EnvironmentCOP passport_env = passports_.top().first.lock();
		if ( passport_env && passport_env->id() == superenv->id() ) {
			throw CREATE_EXCEPTION(EXCN_Env_Passport, "ClientMover being double-assigned a passport for an environment.",
				get_name(), env_ap );
		}
	} else {
		if ( !( passports_.empty() ) ) {
			if ( passports_.top().first.lock()->id() == env->id() ) {
				std::ostringstream ss;
				ss << "ClientMover " << this->get_name() << " already has a passport for " << env->name()
					<< ". This probably means Environment::start got called multiple times." << std::endl;
				throw CREATE_EXCEPTION(EXCN_Env_Passport, ss.str(), get_name(), env_ap );
			} else {
				throw CREATE_EXCEPTION(EXCN_Env_Passport, "ClientMover lacks a superenvironment passport for a superenvironment with which it is registered.",
					get_name(), env_ap );
			}
		}
	}

	passports_.push( std::make_pair( env_ap, pass ) );
}

void ClientMover::pop_passport( Environment const& canceling_env ) {
	// passport_cancel should never be called on an empty passport stack
	if ( passports_.empty() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Environment '"+canceling_env.name()+" trying to pop a passport on "+this->get_name()+"'s empty pass stack." );
	}
	// in fact, we should only ever cancel our own passports. (Should this be an assert?)
	EnvironmentCOP passport_env = passports_.top().first.lock();

	if ( !passport_env ) {
		// If we can't get a lock on our top environment, our top environment has been (or is currently)
		// being deallocated. We'll allow the delete. Optimally, we'd want to check to see that the incoming passport
		// is also inaccessible via owning pointer, but I don't really know how to do that...
	} else if ( passport_env->id() != canceling_env.id() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Environment '"+canceling_env.name()+" trying to pop a passport on "+this->get_name()+" it did not issue." );

	}
	passports_.pop();
	passport_updated();
}

bool ClientMover::state_check( std::string const& method_name, bool test ) const {
	if ( !has_passport() || test ) {
		return true;
	} else {
		std::ostringstream ss;
		ss << "Call to " << this->get_name() << "::" << method_name << " is illegal because it has already "
			<< "issued its claims." << std::endl;
		throw CREATE_EXCEPTION(utility::excn::Exception,  ss.str() );
	}
}


} // environment
} // protocols
