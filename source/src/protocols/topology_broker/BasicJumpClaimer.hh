// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/topology_broker/BasicJumpClaimer.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_topology_broker_BasicJumpClaimer_hh
#define INCLUDED_protocols_topology_broker_BasicJumpClaimer_hh

// Unit headers
#include <protocols/topology_broker/BasicJumpClaimer.fwd.hh>

// Package headers
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <protocols/topology_broker/claims/JumpClaim.hh>

// Project Headers

// Utility headers
#include <utility/vector1.hh>
#include <string>
#include <sstream>

namespace protocols {
namespace topology_broker {

class BasicJumpClaimer : public virtual TopologyClaimer {
	typedef TopologyClaimer Parent;
public:
	virtual ~BasicJumpClaimer() {}

	BasicJumpClaimer() {}

	BasicJumpClaimerOP shared_from_this() { return utility::pointer::dynamic_pointer_cast<BasicJumpClaimer>( TopologyClaimer::shared_from_this() ); }

	TopologyClaimerOP clone() const { return TopologyClaimerOP( new BasicJumpClaimer( *this ) ); }

	std::string type() const { return _static_type_name(); }

	static std::string _static_type_name() { return "BasicJumpClaimer"; }

	virtual void generate_claims( claims::DofClaims& );

	virtual void initialize_dofs( core::pose::Pose&, claims::DofClaims const&, claims::DofClaims&);

	const std::string& end_label() const {
		return end_label_;
	}

	void end_label(const std::string& endLabel) {
		end_label_ = endLabel;
	}

	core::Size end_position() const {
		return end_position_;
	}

	void end_position(core::Size endPosition) {
		end_position_ = endPosition;
	}

	const std::string& start_label() const {
		return start_label_;
	}

	void start_label(const std::string& startLabel) {
		start_label_ = startLabel;
	}

	core::Size start_position() const {
		return start_position_;
	}

	void start_position(core::Size startPosition) {
		start_position_ = startPosition;
	}

protected:
	virtual bool read_tag( std::string, std::istream& );

private:
	std::string start_label_;
	std::string end_label_;
	std::string start_atom_;
	std::string end_atom_;
	core::Size start_position_;
	core::Size end_position_;

}; // class BasicJumpClaimer
} // topology_broker
} // protocols

#endif
