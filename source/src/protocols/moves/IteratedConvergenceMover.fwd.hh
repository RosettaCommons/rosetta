// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/moves/IteratedConvergenceMover.fwd.hh
/// @brief  Forward declaration of the mover class to repeatedly apply a submover until filter convergence is reached
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#ifndef INCLUDED_protocols_moves_IteratedConvergenceMover_fwd_hh
#define INCLUDED_protocols_moves_IteratedConvergenceMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace moves {

class IteratedConvergenceMover;
typedef utility::pointer::shared_ptr< IteratedConvergenceMover > IteratedConvergenceMoverOP;
typedef utility::pointer::shared_ptr< IteratedConvergenceMover const > IteratedConvergenceMoverCOP;

} // moves
} // protocols

#endif
