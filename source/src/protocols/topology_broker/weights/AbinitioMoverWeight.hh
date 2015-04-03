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
///           maintains list of ToplogyClaimers
///           maintains AbinitioMoverWeights -- exclusive or non-exclusively markedup dofs like BackboneClaim, IntraResClaim, JumpClaim
///           generates FoldTree, MoveMap, and collects samplers provided by AbinitioMoverWeights
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_weights_AbinitioMoverWeight_hh
#define INCLUDED_protocols_topology_broker_weights_AbinitioMoverWeight_hh

// Unit Headers

// Package Headers

// Project Headers
#include <core/types.hh>

//#include <core/kinematics/MoveMap.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/pointer/access_ptr.fwd.hh>


//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>


// option key includes


namespace protocols {
namespace topology_broker {
namespace weights {

class AbinitioMoverWeight : public utility::pointer::ReferenceCount {
public:
	virtual ~AbinitioMoverWeight() {};
	virtual core::Real weight( core::Size stageID, core::Real progress /* progress within stage */ ) const = 0;
}; //class AbinitioMoverWeight

// Types
typedef  utility::pointer::shared_ptr< AbinitioMoverWeight >  AbinitioMoverWeightOP;
typedef  utility::pointer::shared_ptr< AbinitioMoverWeight const >  AbinitioMoverWeightCOP;

typedef  utility::pointer::weak_ptr< AbinitioMoverWeight >  AbinitioMoverWeightAP;
typedef  utility::pointer::weak_ptr< AbinitioMoverWeight const >  AbinitioMoverWeightCAP;


}
}
}
#endif
