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


#ifndef INCLUDED_protocols_topology_broker_claims_BBClaim_hh
#define INCLUDED_protocols_topology_broker_claims_BBClaim_hh


// Unit Headers
#include <protocols/topology_broker/claims/BBClaim.fwd.hh>
#include <protocols/topology_broker/claims/DofClaim.hh>


// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/SequenceNumberResolver.hh>

// Project Headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/pointer/ReferenceCount.hh>

//#include <basic/options/option_macros.hh>

//// C++ headers
#include <string>
#include <sstream>


// option key includes


namespace protocols {
namespace topology_broker {
namespace claims {


class BBClaim : public DofClaim {
public:
	BBClaim( TopologyClaimerAP tc, Size pos, ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right)
	{
		local_pos_ = std::make_pair( "DEFAULT", pos );
	}

	BBClaim( TopologyClaimerAP tc, std::pair< std::string, core::Size > local_pos, ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right),
		local_pos_( local_pos )
	{}

	virtual DofClaimOP clone() const { return DofClaimOP( new BBClaim( *this ) ); }

// 	Size get_position() const {
// 		return pos_;
// 	}

	LocalPosition local_position() const {
		return local_pos_;
	}

	core::Size global_position() const {
		TopologyClaimerCOP owner_op( owner() );
		return owner_op->broker().sequence_number_resolver().find_global_pose_number( local_pos_ );
	}

	virtual void toggle( core::kinematics::MoveMap& mm, bool new_setting ) const {
		TopologyClaimerCOP owner_op( owner() );
		core::Size pos = owner_op->broker().sequence_number_resolver().find_global_pose_number( local_pos_ );
		mm.set_bb( pos, new_setting );
	}

	virtual void show(std::ostream& os) const {
		TopologyClaimerCOP owner_op( owner() );
		os << "DofClaim-" << str_type() << " owned by a " << (owner_op ? owner_op->type() : "(Unknown)") << " at ("
			 << local_pos_.first << ", " << local_pos_.second << ")";
	}

//    virtual std::string to_string() const {
//        std::ostringstream str_stream;
//        str_stream << "(BBClaim; owner, " << owner()->type() << "; pos, " << pos_ << ")" ;
//        return str_stream.str();
//    }

	virtual std::string str_type() const {
		return "BB";
	}

protected:
	std::pair< std::string, core::Size > local_pos_;	//first: label, second: position in label

}; //class BBClaim


}
}
}

#endif
