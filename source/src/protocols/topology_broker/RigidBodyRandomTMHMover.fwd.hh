// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rigid/RigidBodyRandomTMHMover.fwd.hh
/// @brief  RigidBodyMover forward declarations header
/// @author

#ifndef INCLUDED_protocols_topology_broker_RigidBodyRandomTMHMover_fwd_hh
#define INCLUDED_protocols_topology_broker_RigidBodyRandomTMHMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace topology_broker{
	
class RigidBodyRandomTMHMover;
typedef utility::pointer::owning_ptr<RigidBodyRandomTMHMover> RigidBodyRandomTMHMoverOP;
typedef utility::pointer::owning_ptr<RigidBodyRandomTMHMover const> RigidBodyRandomTMHMoverCOP;

}
}

#endif