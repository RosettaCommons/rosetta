// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ncbb/ResidueReplacementRebuildMover.fwd.hh
/// @brief A simple method to mutate a pose fairly destructively, and barely recover
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_protocols_ncbb_ResidueReplacementRebuildMover_fwd_hh
#define INCLUDED_protocols_ncbb_ResidueReplacementRebuildMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace ncbb {

class ResidueReplacementRebuildMover;

typedef utility::pointer::shared_ptr< ResidueReplacementRebuildMover > ResidueReplacementRebuildMoverOP;
typedef utility::pointer::shared_ptr< ResidueReplacementRebuildMover const > ResidueReplacementRebuildMoverCOP;

} //protocols
} //ncbb

#endif //INCLUDED_protocols_ncbb_ResidueReplacementRebuildMover_fwd_hh
