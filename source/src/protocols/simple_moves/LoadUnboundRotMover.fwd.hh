// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/LoadUnboundRotMover.fwd.hh
/// @brief
/// @author Steven Lewis smlewi@gmail.com, Brian Weitzner brian.weitzner@gmail.com


#ifndef INCLUDED_protocols_simple_moves_LoadUnboundRotMover_fwd_hh
#define INCLUDED_protocols_simple_moves_LoadUnboundRotMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols{
namespace simple_moves{

//Forwards and OP typedefs
class LoadUnboundRotMover;
typedef utility::pointer::shared_ptr< LoadUnboundRotMover > LoadUnboundRotMoverOP;
typedef utility::pointer::shared_ptr< LoadUnboundRotMover const > LoadUnboundRotMoverCOP;

}//simple_moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_LoadUnboundRotMover_FWD_HH
