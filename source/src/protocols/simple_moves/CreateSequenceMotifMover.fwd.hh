// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/CreateSequenceMotifMover.fwd.hh
/// @brief Create a sequence motif in a region of protein using the SequenceMotifTaskOperation.  Uses psueo-regular expressions to define the motif.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_CreateSequenceMotifMover_fwd_hh
#define INCLUDED_protocols_simple_moves_CreateSequenceMotifMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace simple_moves {

class CreateSequenceMotifMover;

typedef utility::pointer::shared_ptr< CreateSequenceMotifMover > CreateSequenceMotifMoverOP;
typedef utility::pointer::shared_ptr< CreateSequenceMotifMover const > CreateSequenceMotifMoverCOP;

} //protocols
} //simple_moves

#endif //INCLUDED_protocols_simple_moves_CreateSequenceMotifMover_fwd_hh
