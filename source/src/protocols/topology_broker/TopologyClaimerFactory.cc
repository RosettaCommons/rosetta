// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/topology_broker/TopologyClaimerFactory.cc
/// @author Oliver Lange
/// @author Christopher Miles (cmiles@uw.edu)

// Unit Headers
#include <protocols/topology_broker/TopologyClaimerFactory.hh>

// Package Headers
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
#include <protocols/topology_broker/BasicJumpClaimer.hh>

// C/C++ headers
#include <sstream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

static basic::Tracer tr("protocols.topo_broker", basic::t_info);

namespace protocols {
namespace topology_broker {

// Singleton initialization
TopologyClaimerFactory * TopologyClaimerFactory::instance_( 0 );

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex TopologyClaimerFactory::singleton_mutex_;

std::mutex & TopologyClaimerFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
TopologyClaimerFactory const & TopologyClaimerFactory::get_instance()
{
	boost::function< TopologyClaimerFactory * () > creator = boost::bind( &TopologyClaimerFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return *instance_;
}

TopologyClaimerFactory *
TopologyClaimerFactory::create_singleton_instance()
{
	return new TopologyClaimerFactory;
}

// Registers commonly used claimers with the name returned by claimer->type()
TopologyClaimerFactory::TopologyClaimerFactory() {
	add_type(new RigidChunkClaimer());
	add_type(new SequenceClaimer());
	add_type(new FragmentJumpClaimer());
	add_type(new DisulfJumpClaimer());
	add_type(new FragmentClaimer());
	add_type(new ConstraintClaimer());
	add_type(new MembraneTopologyClaimer());
	add_type(new MetalloClaimer());
	add_type(new TemplateJumpClaimer());
	add_type(new CoordConstraintClaimer());
	add_type(new StartStructClaimer());
	add_type(new CutBiasClaimer());
	add_type(new DensityScoringClaimer());
	add_type(new PseudocontactShiftEnergyController());
  add_type(new PseudocontactShiftEnergyController_Ts1());
  add_type(new PseudocontactShiftEnergyController_Ts2());
  add_type(new PseudocontactShiftEnergyController_Ts3());
  add_type(new PseudocontactShiftEnergyController_Ts4());
	add_type(new PcsEnergyController());
	add_type(new FoldandDockClaimer());
	add_type(new FibrilModelingClaimer());
	add_type(new AsymFoldandDockClaimer());
	add_type(new SymmetryClaimer());
	add_type(new BasicJumpClaimer());
}

TopologyClaimerFactory::~TopologyClaimerFactory() {
	delete instance_;
}

void TopologyClaimerFactory::add_type(TopologyClaimerOP claimer) {
	add_type(claimer->type(), claimer);
}

void TopologyClaimerFactory::add_type(const std::string& name, TopologyClaimerOP claimer) {
	claimers_[name] = claimer;
}

TopologyClaimerOP TopologyClaimerFactory::newTopologyClaimer(const std::string& name) const {
	using std::stringstream;

	if (claimers_.find(name) != claimers_.end()) {
		return claimers_[name]->clone();
	} else {
		stringstream ss;
		ss << name
		   << " does not name a known TopologyClaimer -->"
		   << " check spelling or register the type via the add_type() method";
		utility_exit_with_message(ss.str());

		// purely superficial return statement to quiet the compiler
		return NULL;
	}
}

}
}
