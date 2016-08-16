// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/JumpRotamerSidechainMover.fwd.hh
/// @brief
/// @author


#ifndef INCLUDED_protocols_simple_moves_sidechain_moves_JumpRotamerSidechainMover_fwd_hh
#define INCLUDED_protocols_simple_moves_sidechain_moves_JumpRotamerSidechainMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

// Package headers

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

class JumpRotamerSidechainMover;
typedef utility::pointer::shared_ptr< JumpRotamerSidechainMover > JumpRotamerSidechainMoverOP;
typedef utility::pointer::shared_ptr< JumpRotamerSidechainMover const > JumpRotamerSidechainMoverCOP;

} // sidechain_moves
} // simple_moves
} // protocols


#endif  //INCLUDED_protocols_simple_moves_sidechain_moves_JumpRotamerSidechainMover_fwd_HH
