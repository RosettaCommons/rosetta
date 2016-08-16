// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/VectorPoseMover.fwd.hh
/// @brief Used for movers that require multiple poses to act upon
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_PROTOCOLS_MOVES_VectorPoseMover_FWD_HH
#define INCLUDED_PROTOCOLS_MOVES_VectorPoseMover_FWD_HH

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace moves {

class VectorPoseMover;
typedef utility::pointer::shared_ptr<VectorPoseMover> VectorPoseMoverOP;
typedef utility::pointer::shared_ptr<VectorPoseMover const> VectorPoseMoverCOP;

}  // namespace moves
}  // namespace protocols

#endif  // PROTOCOLS_MOVES_CYCLIC_MOVER_FWD_HH_
