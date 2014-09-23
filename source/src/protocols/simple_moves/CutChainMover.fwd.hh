// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/ForceDisulfidesMover.fwd.hh
/// @brief  ForceDisulfidesMover forward declarations header
/// @author Sarel Fleishman

#ifndef INCLUDED_protocols_simple_moves_CutChainMover_fwd_hh
#define INCLUDED_protocols_simple_moves_CutChainMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace simple_moves{

//Forwards and OP typedefs
class CutChainMover;
typedef utility::pointer::shared_ptr< CutChainMover > CutChainMoverOP;
typedef utility::pointer::shared_ptr< CutChainMover const > CutChainMoverCOP;

}//moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_CutChainMover_FWD_HH
