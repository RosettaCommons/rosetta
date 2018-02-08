// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/WriteSSEMover.fwd.hh
/// @brief Writes SSE assignation from DSSP or prediction from PSIPRED as REMARK.
/// @author Jaume Bonet (jaume.bonet@gmail.com)


#ifndef INCLUDED_protocols_simple_filters_WriteSSEMover_fwd_hh
#define INCLUDED_protocols_simple_filters_WriteSSEMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

class WriteSSEMover;
typedef utility::pointer::shared_ptr< WriteSSEMover > WriteSSEMoverOP;
typedef utility::pointer::shared_ptr< WriteSSEMover const > WriteSSEMoverCOP;

}
}

#endif
