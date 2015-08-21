// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/MissingDensityToJumpMover.fwd.hh
/// @brief  Implementation of mover that inserts a jump where there is gap in the pdb. This gap corresponds to missing density.
/// @author TJ Brunette (tjbrunette@gmail.com), May 2011


#ifndef INCLUDED_protocols_simple_moves_MissingDensityToJumpMover_fwd_hh
#define INCLUDED_protocols_simple_moves_MissingDensityToJumpMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

// Package headers

namespace protocols {
namespace simple_moves {

class MissingDensityToJumpMover;
typedef utility::pointer::shared_ptr< MissingDensityToJumpMover > MissingDensityToJumpMoverOP;
typedef utility::pointer::shared_ptr< MissingDensityToJumpMover const > MissingDensityToJumpMoverCOP;

} // moves
} // protocols


#endif
