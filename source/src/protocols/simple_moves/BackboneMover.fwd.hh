// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/BackboneMover.fwd.hh
/// @brief  BackboneMover forward declarations header
/// @author

#ifndef INCLUDED_protocols_simple_moves_BackboneMover_fwd_hh
#define INCLUDED_protocols_simple_moves_BackboneMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

//Forwards and OP typedefs
class BackboneMover;
typedef utility::pointer::shared_ptr< BackboneMover > BackboneMoverOP;
typedef utility::pointer::shared_ptr< BackboneMover const > BackboneMoverCOP;

class SmallMover;
typedef utility::pointer::shared_ptr< SmallMover > SmallMoverOP;
typedef utility::pointer::shared_ptr< SmallMover const > SmallMoverCOP;

class ShearMover;
typedef utility::pointer::shared_ptr< ShearMover > ShearMoverOP;
typedef utility::pointer::shared_ptr< ShearMover const > ShearMoverCOP;

} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_BackboneMover_FWD_HH
