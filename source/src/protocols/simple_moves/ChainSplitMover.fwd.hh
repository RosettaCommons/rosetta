// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/simple_moves/ChainSplitMover.fwd.hh
/// @brief Forward declaration of ChainSplitMover which splits a pose into two chains at a given cutpoint.
/// @author Robert Lindner <rlindner@mpimf-heidelberg.mpg.de>


#ifndef INCLUDED_protocols_simple_moves_ChainSplitMover_FWD_HH
#define INCLUDED_protocols_simple_moves_ChainSplitMover_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace simple_moves {

class ChainSplitMover;

typedef utility::pointer::owning_ptr< ChainSplitMover > ChainSplitMoverOP;
typedef utility::pointer::owning_ptr< ChainSplitMover const > ChainSplitMoverCOP;

}
}


#endif
