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


#ifndef INCLUDED_protocols_topology_broker_claims_LegacyRootClaim_hh
#define INCLUDED_protocols_topology_broker_claims_LegacyRootClaim_hh


// Unit Headers
#include <protocols/topology_broker/claims/LegacyRootClaim.fwd.hh>
#include <protocols/topology_broker/claims/DofClaim.hh>


// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh> //for printing
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>


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

class LegacyRootClaim : public DofClaim {
	//this class could also specify which atoms to use for the jumps --- but I never used this so far.... might be necessary for Zn jumps.
public:
	LegacyRootClaim( TopologyClaimerAP tc, core::Size pos1, ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right ),
		local_position_( "", pos1 ) // An empty label has been hacked to mean no offset in SequenceNumberResolver
	{}

	LegacyRootClaim( TopologyClaimerAP tc, std::pair< std::string, core::Size > local_position, ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right ),
		local_position_ ( local_position )
	{}

	virtual DofClaimOP clone() const { return DofClaimOP( new LegacyRootClaim( *this ) ); }

	std::pair< std::string, core::Size > local_position() const {
		return local_position_;
	}

	virtual void show(std::ostream& os) const {
		os << " with position: " << local_position_.second << " under claim label: " << local_position_.first;
	}

	virtual bool remove() const {
		return false;
	}

	virtual std::string str_type() const {
		return "ROOT";
	}

private:
	// bool permanent_; //true if this cut should still be present after loop-closing
	std::pair< std::string, core::Size > local_position_;
}; //LegacyRootClaim

}
}
}

#endif
