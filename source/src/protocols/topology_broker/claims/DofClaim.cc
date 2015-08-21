// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange

// Unit Headers
//#include <protocols/topology_broker/claims/DofClaim.hh>
//#include <protocols/topology_broker/claims/BBClaim.hh>
//#include <protocols/topology_broker/claims/CutClaim.hh>
//#include <protocols/topology_broker/claims/JumpClaim.hh>
//#include <protocols/topology_broker/claims/LegacyRootClaim.hh>
//#include <protocols/topology_broker/claims/SequenceClaim.hh>
//#include <protocols/topology_broker/claims/SymmetryClaim.hh>

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh> //for printing

// Project Headers

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
//#include <basic/options/option.hh>

//// C++ headers
#include <iostream>

#include <utility/vector1.hh>

// option key includes


static thread_local basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {
namespace claims {

/// @details Auto-generated virtual destructor
DofClaim::~DofClaim() {}

using namespace core;

void DofClaim::show( std::ostream& os ) const {
	os << "owned by, " << owner().lock()->type() << ";";
}

extern std::ostream& operator<<( std::ostream& os, DofClaim const& dof ) {
	dof.show( os );
	return os;
}

extern std::ostream& operator<<( std::ostream& os, DofClaims const& dofs ) {
	for ( DofClaims::const_iterator it = dofs.begin(); it != dofs.end(); ++it ) {
		if ( *it ) {
			os << **it << "\n";
		} else {
			os << "No-Claim\n";
		}
	}
	return os;
}

} //claims
} //topology_broker
} //protocols
