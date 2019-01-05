// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/topology_broker/TopologyClaimerFactory.cc
/// @author Oliver Lange
/// @author Christopher Miles (cmiles@uw.edu)

// Unit Headers
#include <protocols/topology_broker/TopologyClaimerFactory.hh>

// Package Headers
#include <protocols/topology_broker/TopologyBroker.fwd.hh>
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>
#include <protocols/topology_broker/RigidChunkClaimer.hh>
#include <protocols/topology_broker/ConstraintClaimer.hh>
#include <protocols/topology_broker/FragmentClaimer.hh>
#include <protocols/topology_broker/FragmentJumpClaimer.hh>
#include <protocols/topology_broker/DisulfJumpClaimer.hh>
#include <protocols/topology_broker/SequenceClaimer.hh>
#include <protocols/topology_broker/SymmetryClaimer.hh>
#include <protocols/topology_broker/MetalloClaimer.hh>
#include <protocols/topology_broker/MembraneTopologyClaimer.hh>
#include <protocols/topology_broker/TemplateJumpClaimer.hh>
#include <protocols/topology_broker/CoordConstraintClaimer.hh>
#include <protocols/topology_broker/StartStructClaimer.hh>
#include <protocols/topology_broker/CutBiasClaimer.hh>
#include <protocols/topology_broker/DensityScoringClaimer.hh>
#include <protocols/topology_broker/FoldandDockClaimer.hh>
#include <protocols/topology_broker/AsymFoldandDockClaimer.hh>
#include <protocols/topology_broker/PseudocontactShiftEnergyController.hh>
#include <protocols/topology_broker/PseudocontactShiftEnergyController_Ts1.hh>
#include <protocols/topology_broker/PseudocontactShiftEnergyController_Ts2.hh>
#include <protocols/topology_broker/PseudocontactShiftEnergyController_Ts3.hh>
#include <protocols/topology_broker/PseudocontactShiftEnergyController_Ts4.hh>
#include <protocols/topology_broker/PcsEnergyController.hh>
#include <protocols/topology_broker/FibrilModelingClaimer.hh>
#include <protocols/topology_broker/TMHTopologySamplerClaimer.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <protocols/topology_broker/BasicJumpClaimer.hh>

// C/C++ headers
#include <sstream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

static basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

// Registers commonly used claimers with the name returned by claimer->type()
TopologyClaimerFactory::TopologyClaimerFactory() {
	add_type(utility::pointer::make_shared< RigidChunkClaimer >());
	add_type(utility::pointer::make_shared< SequenceClaimer >());
	add_type(utility::pointer::make_shared< FragmentJumpClaimer >());
	add_type(utility::pointer::make_shared< DisulfJumpClaimer >());
	add_type(utility::pointer::make_shared< FragmentClaimer >());
	add_type(utility::pointer::make_shared< ConstraintClaimer >());
	add_type(utility::pointer::make_shared< MembraneTopologyClaimer >());
	add_type(utility::pointer::make_shared< MetalloClaimer >());
	add_type(utility::pointer::make_shared< TemplateJumpClaimer >());
	add_type(utility::pointer::make_shared< CoordConstraintClaimer >());
	add_type(utility::pointer::make_shared< StartStructClaimer >());
	add_type(utility::pointer::make_shared< CutBiasClaimer >());
	add_type(utility::pointer::make_shared< DensityScoringClaimer >());
	add_type(utility::pointer::make_shared< PseudocontactShiftEnergyController >());
	add_type(utility::pointer::make_shared< PseudocontactShiftEnergyController_Ts1 >());
	add_type(utility::pointer::make_shared< PseudocontactShiftEnergyController_Ts2 >());
	add_type(utility::pointer::make_shared< PseudocontactShiftEnergyController_Ts3 >());
	add_type(utility::pointer::make_shared< PseudocontactShiftEnergyController_Ts4 >());
	add_type(utility::pointer::make_shared< PcsEnergyController >());
	add_type(utility::pointer::make_shared< FoldandDockClaimer >());
	add_type(utility::pointer::make_shared< FibrilModelingClaimer >());
	add_type(utility::pointer::make_shared< AsymFoldandDockClaimer >());
	add_type(utility::pointer::make_shared< TMHTopologySamplerClaimer >());
	add_type(utility::pointer::make_shared< SymmetryClaimer >());
	add_type(utility::pointer::make_shared< BasicJumpClaimer >());
}

TopologyClaimerFactory::~TopologyClaimerFactory() = default;

void TopologyClaimerFactory::add_type(TopologyClaimerOP claimer) {
	add_type(claimer->type(), claimer);
}

void TopologyClaimerFactory::add_type(const std::string& name, TopologyClaimerOP claimer) {
	claimers_[name] = claimer;
}

TopologyClaimerOP TopologyClaimerFactory::newTopologyClaimer(const std::string& name) const {
	using std::stringstream;

	if ( claimers_.find(name) != claimers_.end() ) {
		return claimers_[name]->clone();
	} else {
		stringstream ss;
		ss << name
			<< " does not name a known TopologyClaimer -->"
			<< " check spelling or register the type via the add_type() method";
		utility_exit_with_message(ss.str());

		// purely superficial return statement to quiet the compiler
		return nullptr;
	}
}

}
}
