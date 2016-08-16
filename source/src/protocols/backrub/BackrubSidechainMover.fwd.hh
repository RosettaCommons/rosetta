// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/backrub/BackrubSidechainMover.fwd.hh
/// @brief  BackrubSidechainMover forward declarations header
/// @author

#ifndef INCLUDED_protocols_backrub_BackrubSidechainMover_fwd_hh
#define INCLUDED_protocols_backrub_BackrubSidechainMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace backrub {

// Forwards and OP typedefs
class BackrubSidechainMover;
typedef utility::pointer::shared_ptr< BackrubSidechainMover > BackrubSidechainMoverOP;
typedef utility::pointer::shared_ptr< BackrubSidechainMover const > BackrubSidechainMoverCOP;

} //moves
} //protocols

#endif //INCLUDED_protocols_backrub_BackrubSidechainMover_FWD_HH
