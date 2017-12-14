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
	add_type(TopologyClaimerOP( new RigidChunkClaimer() ));
	add_type(TopologyClaimerOP( new SequenceClaimer() ));
	add_type(TopologyClaimerOP( new FragmentJumpClaimer() ));
	add_type(TopologyClaimerOP( new DisulfJumpClaimer() ));
	add_type(TopologyClaimerOP( new FragmentClaimer() ));
	add_type(TopologyClaimerOP( new ConstraintClaimer() ));
	add_type(TopologyClaimerOP( new MembraneTopologyClaimer() ));
	add_type(TopologyClaimerOP( new MetalloClaimer() ));
	add_type(TopologyClaimerOP( new TemplateJumpClaimer() ));
	add_type(TopologyClaimerOP( new CoordConstraintClaimer() ));
	add_type(TopologyClaimerOP( new StartStructClaimer() ));
	add_type(TopologyClaimerOP( new CutBiasClaimer() ));
	add_type(TopologyClaimerOP( new DensityScoringClaimer() ));
	add_type(TopologyClaimerOP( new PseudocontactShiftEnergyController() ));
	add_type(TopologyClaimerOP( new PseudocontactShiftEnergyController_Ts1() ));
	add_type(TopologyClaimerOP( new PseudocontactShiftEnergyController_Ts2() ));
	add_type(TopologyClaimerOP( new PseudocontactShiftEnergyController_Ts3() ));
	add_type(TopologyClaimerOP( new PseudocontactShiftEnergyController_Ts4() ));
	add_type(TopologyClaimerOP( new PcsEnergyController() ));
	add_type(TopologyClaimerOP( new FoldandDockClaimer() ));
	add_type(TopologyClaimerOP( new FibrilModelingClaimer() ));
	add_type(TopologyClaimerOP( new AsymFoldandDockClaimer() ));
	add_type(TopologyClaimerOP( new TMHTopologySamplerClaimer() ));
	add_type(TopologyClaimerOP( new SymmetryClaimer() ));
	add_type(TopologyClaimerOP( new BasicJumpClaimer() ));
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
