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
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_LoopFragmentClaimer_hh
#define INCLUDED_protocols_topology_broker_LoopFragmentClaimer_hh


// Unit Headers
#include <protocols/topology_broker/LoopFragmentClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/FragmentClaimer.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.hh>
//#include <core/fragment/LoopFragSet.hh>

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
//#include <basic/Tracer.hh>
//#include <basic/options/option.hh>

#include <utility/pointer/ReferenceCount.hh>

//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>
#include <string>

#include <utility/vector1.hh>


// option key includes


namespace protocols {
namespace topology_broker {

class LoopFragmentClaimer : public FragmentClaimer {
public:
	LoopFragmentClaimer( core::fragment::FragSetOP );
	LoopFragmentClaimer( core::fragment::FragSetOP, std::string label );


	/// @brief type() is specifying the output name of the TopologyClaimer
	std::string type() const override {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "LoopFragmentClaimer";
	}

	moves::MoverOP get_mover( core::pose::Pose const& /*pose*/ ) const override {
		return nullptr; //override FragmentClaimer functionality
	}

	core::fragment::FragSetCOP loop_frags( core::kinematics::MoveMap& movemap) const override;
	// { virtual void loop_frags( core::fragment::FragSet& ) const;

}; //class LoopFragmentClaimer


}
}

#endif
