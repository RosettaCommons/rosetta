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
///           maintains AbinitioMoverWeights -- exclusive or non-exclusively markedup dofs like BackboneClaim, IntraResClaim, JumpClaim
///           generates FoldTree, MoveMap, and collects samplers provided by AbinitioMoverWeights
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_weights_ConstAbinitioMoverWeight_hh
#define INCLUDED_protocols_topology_broker_weights_ConstAbinitioMoverWeight_hh

// Unit Headers
#include <protocols/topology_broker/weights/AbinitioMoverWeight.hh>

// Package Headers

// Project Headers
#include <core/types.hh>

//#include <core/kinematics/MoveMap.fwd.hh>

// ObjexxFCL Headers

// Utility headers

#include <utility/pointer/ReferenceCount.hh>

//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>


// option key includes


namespace protocols {
namespace topology_broker {
namespace weights {


class ConstAbinitioMoverWeight : public AbinitioMoverWeight {
public:
	ConstAbinitioMoverWeight( core::Real weight ) : weight_( weight ) {};
	virtual core::Real weight( core::Size /* stageID */, core::Real /* progress */ /* progress within stage */ )  const {
		return weight_;
	}
protected:
	core::Real base_weight() { return weight_; }
private:
	core::Real weight_;
}; //class ConstAbinitioMoverWeight

// Types
typedef  utility::pointer::shared_ptr< ConstAbinitioMoverWeight >  ConstAbinitioMoverWeightOP;
typedef  utility::pointer::shared_ptr< ConstAbinitioMoverWeight const >  ConstAbinitioMoverWeightCOP;

typedef  utility::pointer::weak_ptr< ConstAbinitioMoverWeight >  ConstAbinitioMoverWeightAP;
typedef  utility::pointer::weak_ptr< ConstAbinitioMoverWeight const >  ConstAbinitioMoverWeightCAP;

}
}
}

#endif
