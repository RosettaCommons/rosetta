// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/ClientMover.hh
/// @brief An environment that automatically distributes rights to shared degrees of freedom (e.g. fold tree)
///
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_ClientMover_hh
#define INCLUDED_protocols_environment_ClientMover_hh

// Unit Headers
#include <protocols/environment/ClientMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <protocols/environment/Environment.hh>

#include <protocols/environment/claims/EnvClaim.fwd.hh>
#include <protocols/environment/claims/VirtResClaim.fwd.hh>

#include <core/environment/DofPassport.fwd.hh>

// Project headers
#include <core/conformation/Conformation.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/WriteableCacheableMap.fwd.hh>

// C++ Headers
#include <stack>

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class ClientMover : public protocols::moves::Mover {
	// Uses push and pop_passport.
	friend class Environment;
	// Uses passport_updated during initialization
	friend class EnvClaimBroker;

	typedef core::environment::DofPassportCOP DofPassportCOP;

	typedef std::stack< std::pair< EnvironmentCAP, DofPassportCOP > > PassportStack;

public:
	ClientMover();

	ClientMover( ClientMover const& );

	virtual ~ClientMover();

	/// @brief   Returns a list of claims for this mover.
	/// @details The pose passed as an argument is used for reference informational
	///          purposes only (for example, you want to know the sequence) or access
	///          the pose cache. Any changes to the conformation will get overwritten!
	virtual claims::EnvClaims yield_claims( core::pose::Pose const&,
		basic::datacache::WriteableCacheableMapOP ) = 0;

	/// @brief this method is called by the broking system in response to a
	///        successful initialization claim by this mover. The passport
	///        will prevent any unauthorized changes.
	virtual void initialize( Pose& conf );

	/// @brief this method is used to make sure any movers contained in this
	///        mover (that want to make claims) also get registered with the
	///        environment.
	virtual void yield_submovers( std::set< ClientMoverOP >& ) const {};

protected:

	/// @brief convience method that copies all changes from an unprotected pose to a protected pose.
	void sandboxed_copy( core::pose::Pose const& sandbox_pose,
		core::pose::Pose& true_pose ) const;

	/// @brief hook that is called each time the passport status is updated
	///        (presumably because there's a new subenvironment, or the sub-
	///        environment has expired).
	virtual void passport_updated() {}

	/// @brief hook that provides information about the final result of broking
	///        at the very end, though this BrokerResult construct.
	virtual void broking_finished( EnvClaimBroker::BrokerResult const& ) {}

	/// @brief return a pointer the active environment
	EnvironmentCAP active_environment() const;

	/// @brief read-only access to the top passport in the passport stack.
	core::environment::DofPassportCOP passport() const;

	bool has_passport() const {
		return !passports_.empty();
	}

	/// @brief convienence method for failing a state-set call after brokering.
	/// @details once brokering has been completed, most state-set calls
	///          (e.g. set_mover or the like) don't make sense anymore. This
	///          method is used to convienently throw an explanatory exception.
	bool state_check( std::string const& method_name, bool test ) const;

private:

	/// @brief called by environments to push a new passport on the passport stack
	//friend void Environment::assign_passport( ClientMoverOP, DofPassportCOP );

	/// @brief called by environments when the environment expires to pop the environment.
	//friend void Environment::cancel_passport( ClientMoverOP );

	void push_passport( EnvironmentCAP, DofPassportCOP );

	void pop_passport( Environment const& );

	PassportStack passports_;

}; // end ClientMover base class

} // moves
} // protocols

#endif //INCLUDED_protocols_environment_ClientMover_HH
