// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/RandomTorsionMover.fwd.hh
/// @brief  RandomTorsionMover forward declarations header
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

#ifndef INCLUDED_protocols_simple_moves_RandomTorsionMover_FWD_HH
#define INCLUDED_protocols_simple_moves_RandomTorsionMover_FWD_HH

// utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

class RandomTorsionMover;

typedef utility::pointer::shared_ptr< RandomTorsionMover > RandomTorsionMoverOP;
typedef utility::pointer::shared_ptr< RandomTorsionMover const > RandomTorsionMoverCOP;

} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_RandomTorsionMover_FWD_HH
