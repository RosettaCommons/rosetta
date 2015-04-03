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


#ifndef INCLUDED_protocols_topology_broker_claims_CutClaim_hh
#define INCLUDED_protocols_topology_broker_claims_CutClaim_hh


// Unit Headers
#include <protocols/topology_broker/claims/CutClaim.fwd.hh>
#include <protocols/topology_broker/claims/DofClaim.hh>


// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh> //for printing
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
//#include <basic/Tracer.hh>
//#include <basic/options/option.hh>
#include <utility/exit.hh>
#include <utility/pointer/ReferenceCount.hh>

//#include <basic/options/option_macros.hh>

#include <utility/vector1.hh>


//// C++ headers
#include <utility>
#include <string>
#include <sstream>


// option key includes


namespace protocols {
namespace topology_broker {
namespace claims {

class CutClaim : public DofClaim {
	//this class could also specify which atoms to use for the jumps --- but I never used this so far.... might be necessary for Zn jumps.
public:
	CutClaim( TopologyClaimerAP tc, std::pair< std::string, core::Size > position,
						ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right ),
		position_( position )
	{}

	virtual DofClaimOP clone() const { return DofClaimOP( new CutClaim( *this ) ); }

	std::pair< std::string, core::Size > get_position() const {
		return position_;
	}

	virtual void show(std::ostream& os) const {
		TopologyClaimerCOP owner_op( owner() );
		os << "DofClaim-" << str_type() << " owned by a " << (owner_op ? owner_op->type() : "(Unknown)") << "  at pos";
		os << " ('" << position_.first << "'," << position_.second << ")";
	}

	virtual bool remove() const {
		return false; //so far all CutClaims will be physical --> permanent cuts ... !permanent_;
	}

	virtual std::string str_type() const {
		return "CUT";
	}

//    virtual std::string to_string() const {
//        std::ostringstream str_stream;
//        str_stream << "(CutClaim; owner, " << owner()->type() << "; pos, " << pos1_ << ")" ;
//        return str_stream.str();
//    }

private:
	//	bool permanent_; //true if this cut should still be present after loop-closing
	std::pair< std::string, core::Size > position_;
}; //CutClaim


}
}
}

#endif


