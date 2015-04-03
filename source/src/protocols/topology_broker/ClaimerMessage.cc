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
#include <protocols/topology_broker/ClaimerMessage.hh>
#include <protocols/topology_broker/TopologyClaimer.hh>
#include <protocols/topology_broker/claims/DofClaim.hh>

#include <iostream>

#include <utility/vector1.hh>


namespace protocols {
namespace topology_broker {

ClaimerMessage::ClaimerMessage() : label_( "NoAddress" ) {}
ClaimerMessage::ClaimerMessage( std::string label ) : label_( label ) {}

ClaimerMessage::~ClaimerMessage() {}


extern std::ostream& operator << ( std::ostream& os, ClaimerMessage const& cm ) {
	os << "ClaimerMessage type: " << cm.type() << " received by:\n ";
	for ( ClaimerMessage::TopologyClaimerCAPs::const_iterator it = cm.received_by_.begin(); it != cm.received_by_.end(); ++it ) {
		TopologyClaimerCOP claim(*it);
		os << "      " << claim->label() << " of type " << claim->type() << "\n";
	}
	cm.show( os );
	return os;
}

void ClaimerMessage::received_by_me( TopologyClaimerCAP me ) {
	received_by_.push_back( me );
}

core::Size
ClaimerMessage::nr_recepients() const {
	return received_by_.size();
}

bool ClaimerMessage::received() const {
	return nr_recepients() > 0;
}

SuggestValueMessage::SuggestValueMessage( std::string label ) : ClaimerMessage( label ) {}
SuggestValueMessage::~SuggestValueMessage() {}


} //topology_broker
} //protocols
