// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/sidechain_moves/SetChiMover.fwd.hh
/// @brief  A mover to change one chi angle
/// @author Noah Ollikanen

#ifndef INCLUDED_protocols_simple_moves_sidechain_moves_SetChiMover_fwd_hh
#define INCLUDED_protocols_simple_moves_sidechain_moves_SetChiMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

// Package headers

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

class SetChiMover;
typedef utility::pointer::shared_ptr< SetChiMover > SetChiMoverOP;
typedef utility::pointer::shared_ptr< SetChiMover const > SetChiMoverCOP;

} // sidechain_moves
} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_sidechain_moves_SetChiMover_FWD_HH
