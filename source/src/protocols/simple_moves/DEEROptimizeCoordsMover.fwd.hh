// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/simple_moves/DEEROptimizeCoordsMover.hh
/// @brief
/// @author Diego del Alamo

#ifndef INCLUDED_protocols_simple_moves_DEEROptimizeCoordsMover_fwd_hh
#define INCLUDED_protocols_simple_moves_DEEROptimizeCoordsMover_fwd_hh

#include <utility/pointer/std/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

class DEEROptimizeCoordsMover;

typedef utility::pointer::shared_ptr< DEEROptimizeCoordsMover > DEEROptimizeCoordsMoverOP;
typedef utility::pointer::shared_ptr< DEEROptimizeCoordsMover const > DEEROptimizeCoordsMoverCOP;

}
}

#endif //INCLUDED_protocols_simple_moves_MSDMover_fwd_hh
