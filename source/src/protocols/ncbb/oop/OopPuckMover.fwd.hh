// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ncbb/oop/OopPuckMover.fwd.hh
/// @brief  OopPuckMover forward declarations header
/// @author Kevin Drew, kdrew@nyu.edu
#ifndef INCLUDED_protocols_simple_moves_oop_OopPuckMover_fwd_hh
#define INCLUDED_protocols_simple_moves_oop_OopPuckMover_fwd_hh
// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {
namespace oop {


//Forwards and OP typedefs
/*
class OopPuckMover;
typedef utility::pointer::owning_ptr< OopPuckMover > OopPuckMoverOP;
typedef utility::pointer::owning_ptr< OopPuckMover const > OopPuckMoverCOP;
*/

class OopPuckPlusMover;
typedef utility::pointer::shared_ptr< OopPuckPlusMover > OopPuckPlusMoverOP;
typedef utility::pointer::shared_ptr< OopPuckPlusMover const > OopPuckPlusMoverCOP;

class OopPuckMinusMover;
typedef utility::pointer::shared_ptr< OopPuckMinusMover > OopPuckMinusMoverOP;
typedef utility::pointer::shared_ptr< OopPuckMinusMover const > OopPuckMinusMoverCOP;

class OopDPuckPlusMover;
typedef utility::pointer::shared_ptr< OopDPuckPlusMover > OopDPuckPlusMoverOP;
typedef utility::pointer::shared_ptr< OopDPuckPlusMover const > OopDPuckPlusMoverCOP;

class OopDPuckMinusMover;
typedef utility::pointer::shared_ptr< OopDPuckMinusMover > OopDPuckMinusMoverOP;
typedef utility::pointer::shared_ptr< OopDPuckMinusMover const > OopDPuckMinusMoverCOP;

}//oop
}//simple_moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_oop_OopPuckMover_fwd_hh
