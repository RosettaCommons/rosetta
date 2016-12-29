// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/PyMOLMover.hh
/// @brief  Send infromation to PyMOL
/// @author Sergey Lyskov

#ifndef INCLUDED_protocols_moves_PyMOLMover_fwd_hh
#define INCLUDED_protocols_moves_PyMOLMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace moves {

class PyMOLMover;
typedef utility::pointer::shared_ptr< PyMOLMover > PyMOLMoverOP;
typedef utility::pointer::shared_ptr< PyMOLMover const > PyMOLMoverCOP;

class PyMOLObserver;
typedef utility::pointer::shared_ptr< PyMOLObserver > PyMOLObserverOP;
typedef utility::pointer::shared_ptr< PyMOLObserver const > PyMOLObserverCOP;

} // moves
} // protocols

#endif // INCLUDED_protocols_moves_PyMOLMover_FWD_HH
