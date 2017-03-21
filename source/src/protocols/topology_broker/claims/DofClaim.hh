// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DofClaim
/// @brief
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_claims_DofClaim_hh
#define INCLUDED_protocols_topology_broker_claims_DofClaim_hh


// Unit Headers
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>

// Package Headers
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

//// C++ headers
//#include <fstream>
#include <string>

#include <utility/vector1.hh>


// option key includes


namespace protocols {
namespace topology_broker {
namespace claims {
/// A better DofClaims class would provide some extracting functions:
/// by owner
/// by type

typedef std::pair< std::string, core::Size> LocalPosition;

class DofClaim : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~DofClaim();
	typedef core::Size Size;

	enum ClaimRight {
		NEED_TO_KNOW = 1,
		CAN_INIT,
		INIT,
		EXCLUSIVE, //MUST INIT
		REJECTED
	};

	DofClaim( TopologyClaimerAP tc, ClaimRight right ) :
		claim_source_( tc ),
		right_( right ),
		approved_( false )
	{};

	virtual DofClaimOP clone() const = 0;
	//virtual?

	ClaimRight right() const { return right_; };

	TopologyClaimerCAP owner() const { return claim_source_; }
	TopologyClaimerAP owner() { return claim_source_; }

	virtual void toggle( core::kinematics::MoveMap&, bool /*new_setting*/ ) const {};

	bool exclusive() const {
		return right() == DofClaim::EXCLUSIVE;
	}

	core::Size last_residue() const {
		utility_exit_with_message("DofClaim::last_residue() is currently unimplemented.");
		return 0;
	}
	//virtual std::string to_string() const = 0;
	virtual std::string str_type() const = 0;
	virtual void show( std::ostream& os ) const;

	bool approved() const {
		return approved_;
	}

	void set_approved() { //this should only be called by the TopologyBroker ... make sure somehow? friend ?
		approved_ = true;
	}
private:
	TopologyClaimerAP claim_source_; //NEVER Make this OP --- circularity in smart-pointers   (wanted this reference but there was some kind of problem ... what was it ?
	ClaimRight right_;
	bool approved_; //keep track of this ?
}; //class DofClaim

extern std::ostream& operator<<( std::ostream& os, DofClaim const& );
extern std::ostream& operator<<( std::ostream& os, DofClaims const& );

}
}
}

#endif
