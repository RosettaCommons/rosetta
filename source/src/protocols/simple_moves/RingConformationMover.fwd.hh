// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/simple_moves/RingConformationMover.fwd.hh
/// @brief   Forward declarations for RingConformationMover.
/// @author  Labonte <JWLabonte@jhu.edu>

#ifndef INCLUDED_protocols_simple_moves_RingConformationMover_FWD_HH
#define INCLUDED_protocols_simple_moves_RingConformationMover_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

/// @brief  A Mover class for switching among various ring conformations.
class RingConformationMover;

typedef utility::pointer::shared_ptr< RingConformationMover > RingConformationMoverOP;
typedef utility::pointer::shared_ptr< RingConformationMover const > RingConformationMoverCOP;

}  // namespace simple_moves
}  // namespace protocols

#endif  // INCLUDED_protocols_simple_moves_RingConformationMover_FWD_HH
