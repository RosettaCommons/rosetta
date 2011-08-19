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
/// @detailed responsibilities:
/// @author Oliver Lange

// Unit Headers
#include <protocols/topology_broker/DofClaim.hh>

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh> //for printing

// Project Headers

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
#include <basic/Tracer.hh>
//#include <basic/options/option.hh>

//// C++ headers
#include <iostream>

//#include <fstream>

// option key includes


static basic::Tracer tr("protocols.topo_broker",basic::t_info);

namespace protocols {
namespace topology_broker {

using namespace core;

void DofClaim::show( std::ostream& os ) const {
	os << "DofClaim-" << str_type() << " owned by a " << owner()->type() << "  at pos";
	for ( Size i = 1; i<=size(); ++i ) { os << " " << pos(i);}
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

} //topology_broker
} //protocols
