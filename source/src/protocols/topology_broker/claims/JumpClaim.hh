// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_claims_JumpClaim_hh
#define INCLUDED_protocols_topology_broker_claims_JumpClaim_hh


// Unit Headers
#include <protocols/topology_broker/claims/JumpClaim.fwd.hh>
#include <protocols/topology_broker/claims/DofClaim.hh>



// Package Headers
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/TopologyClaimer.hh>
#include <protocols/topology_broker/SequenceNumberResolver.hh>

// Project Headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/exit.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>


//// C++ headers
#include <string>
#include <sstream>


// option key includes


namespace protocols {
namespace topology_broker {
namespace claims {

class JumpClaim : public DofClaim {
	//this class could also specify which atoms to use for the jumps --- but I never used this so far.... might be necessary for Zn jumps.
public:

	JumpClaim( TopologyClaimerAP tc,
						 core::Size pos1,
						 core::Size pos2,
						 std::string atom1,
						 std::string atom2,
						 ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right ),
		permanent_( false ),
		atom1_( atom1 ),
		atom2_( atom2 )
	{
		local_pos1_ = std::make_pair( "DEFAULT", pos1 );
		local_pos2_ = std::make_pair( "DEFAULT", pos2 );
	}

	JumpClaim( TopologyClaimerAP tc,
						 core::Size pos1,
						 core::Size pos2,
						 ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right ),
		permanent_( false ),
		atom1_( "" ),
		atom2_( "" )
	{
		local_pos1_ = std::make_pair( "DEFAULT", pos1 );
		local_pos2_ = std::make_pair( "DEFAULT", pos2 );
	}

	JumpClaim( TopologyClaimerAP tc,
						 LocalPosition pos1,
						 LocalPosition pos2,
						 ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right ),
		permanent_( false ),
		local_pos1_( pos1 ),
		local_pos2_( pos2 ),
		atom1_( "" ),
		atom2_( "" )
	{}

	JumpClaim( TopologyClaimerAP tc,
						 LocalPosition pos1,
						 LocalPosition pos2,
						 std::string atom1,
						 std::string atom2,
						 ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right ),
		permanent_( false ),
		local_pos1_( pos1 ),
		local_pos2_( pos2 ),
		atom1_( atom1 ),
		atom2_( atom2 )
	{}

	virtual DofClaimOP clone() const { return DofClaimOP( new JumpClaim( *this ) ); }

	LocalPosition const& local_pos1() const {
		return local_pos1_;
	}

	LocalPosition const& local_pos2() const {
		return local_pos2_;
	}

	virtual void show(std::ostream& os) const {
		TopologyClaimerCOP owner_op( owner() );
		os << "DofClaim-" << str_type() << " owned by a " << (owner_op ? owner_op->type() : "(Unknown)") << " from ("
			 << local_pos1_.first << ", " << local_pos1_.second << ") to ("
			 << local_pos2_.first   << ", " << local_pos2_.second <<").";
	}

	core::Size global_pos1() const {
		TopologyClaimerCOP owner_op( owner() );
		return owner_op->broker().sequence_number_resolver().find_global_pose_number( local_pos1_);
	}

	core::Size global_pos2() const {
		TopologyClaimerCOP owner_op( owner() );
		return owner_op->broker().sequence_number_resolver().find_global_pose_number( local_pos2_ );
	}

	virtual void toggle( core::kinematics::MoveMap& mm, bool new_setting ) const {
		TopologyClaimerCOP owner_op( owner() );
		core::Size pos1 = owner_op->broker().sequence_number_resolver().find_global_pose_number( local_pos1_ );
		core::Size pos2 = owner_op->broker().sequence_number_resolver().find_global_pose_number( local_pos2_ );
		mm.set_jump( pos1, pos2, new_setting );
	}


	virtual bool remove() const {
		return !permanent_;
	}

	virtual std::string str_type() const {
		return "JUMP";
	}

	std::string const& jump_atom1() const {
		return atom1_;
	}

	std::string const& jump_atom2() const {
		return atom2_;
	}

	void jump_atom1( std::string const& str ) {
		atom1_ = str;
	}

	void jump_atom2( std::string const& str ) {
		atom2_ = str;
	}

	std::string const& jump_atom( Size i ) const {
		runtime_assert( i <= 2 && i > 0 );
		if ( i == 1 ) return jump_atom1();
		else return jump_atom2();
	}

	void set_jump_atom( core::Size i, std::string const& str ) {
		runtime_assert( i <= 2 && i > 0 );
		if ( i == 1 ) jump_atom1( str );
		else jump_atom2( str );
	}
  
private:
	bool permanent_; //true if this jump should still be present after loop-closing
// 	Size pos1_;
// 	Size pos2_;
	LocalPosition local_pos1_;
	LocalPosition local_pos2_;
	std::string atom1_;
	std::string atom2_;
}; //JumpClaim

}
}
}

#endif

