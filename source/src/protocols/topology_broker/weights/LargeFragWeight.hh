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
///           maintains SmallFragWeights -- exclusive or non-exclusively markedup dofs like BackboneClaim, IntraResClaim, JumpClaim
///           generates FoldTree, MoveMap, and collects samplers provided by LargeFragWeights
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_weights_LargeFragWeight_hh
#define INCLUDED_protocols_topology_broker_weights_LargeFragWeight_hh

// Unit Headers
#include <protocols/topology_broker/weights/AbinitioMoverWeight.hh>

// Package Headers
#include <protocols/abinitio/FragmentSampler.fwd.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/broker.OptionKeys.gen.hh>

// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>


namespace protocols {
namespace topology_broker {
namespace weights {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace basic::options;
using namespace basic::options::OptionKeys;

class LargeFragWeight : public AbinitioMoverWeight {
public:
	LargeFragWeight( core::Real weight = 1.0 ) : weight_( weight ) {};
	virtual core::Real weight( core::Size stageID, core::Real /*progress*/ /* progress within stage */ ) const
	{
		if ( option[basic::options::OptionKeys::broker::large_frag_mover_stage1_weight].user() && stageID == 1 ) {
			return option[basic::options::OptionKeys::broker::large_frag_mover_stage1_weight].value();
		}

		if ( stageID < abinitio::STAGE_4 ) return weight_;
		return 0.0;
	};
private:
	core::Real weight_;
}; //class LargeFragWeight

// Types
typedef  utility::pointer::shared_ptr< LargeFragWeight >  LargeFragWeightOP;
typedef  utility::pointer::shared_ptr< LargeFragWeight const >  LargeFragWeightCOP;

typedef  utility::pointer::weak_ptr< LargeFragWeight >  LargeFragWeightAP;
typedef  utility::pointer::weak_ptr< LargeFragWeight const >  LargeFragWeightCAP;

}
}
}

#endif
