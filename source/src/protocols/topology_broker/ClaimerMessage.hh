// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
///           maintains list of ToplogyClaimers
///           maintains ClaimerMessages -- exclusive or non-exclusively markedup dofs like BackboneClaim, IntraResClaim, JumpClaim
///           generates FoldTree, MoveMap, and collects samplers provided by ClaimerMessages
/// @author Oliver Lange
#ifndef INCLUDED_protocols_topology_broker_ClaimerMessage_hh
#define INCLUDED_protocols_topology_broker_ClaimerMessage_hh


// Unit Headers
#include <protocols/topology_broker/ClaimerMessage.fwd.hh>
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>

// Package Headers

// Project Headers
#include <core/types.hh>

//include <core/kinematics/MoveMap.fwd.hh>

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
//#include <basic/Tracer.hh>
//#include <basic/options/option.hh>
//#include <utility/fix_boinc_read.hh>

//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>
#include <string>
#include <typeinfo>

#include <utility/vector1_bool.hh>


// option key includes


namespace protocols {
namespace topology_broker {


class ClaimerMessage { //: public utility::pointer::ReferenceCount {
	typedef utility::vector1< TopologyClaimerCAP > TopologyClaimerCAPs;
public:
	ClaimerMessage();
	ClaimerMessage( std::string const & label );
	virtual ~ClaimerMessage();

	/// @brief name of Claimer
	virtual std::string type() const {
		return typeid( *this ).name();
	}

	std::string const& label() const {
		return label_;
	}

	void set_label( std::string const& str ) {
		label_ = str;
	}

	bool matches( std::string const& label ) {
		if ( label_ == "ALL" ) return true;
		return label == label_;
	}

	virtual void show( std::ostream& ) const {};

	void received_by_me( TopologyClaimerCAP me );

	core::Size nr_recepients() const;

	bool received() const;

	friend std::ostream& operator << ( std::ostream& os, ClaimerMessage const& cm );

private:
	/// @brief a user defined string, can be used to send messages from claimer to claimer
	std::string label_;

	TopologyClaimerCAPs received_by_;

}; //class ClaimerMessage


class SuggestValueMessage : public ClaimerMessage {
public:
	SuggestValueMessage( std::string label );
	~SuggestValueMessage() override;
	claims::DofClaimOP some_claim_; //e.g., a ROOT Claim
};


}
}

#endif
