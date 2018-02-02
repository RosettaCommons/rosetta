// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/AlignByResidueSelectorMover.fwd.hh
/// @brief Aligns two poses through the selected residues
/// @author Jaume Bonet (jaume.bonet@gmail.com)


#ifndef INCLUDED_protocols_fold_from_loops_movers_AlignByResidueSelectorMover_fwd_hh
#define INCLUDED_protocols_fold_from_loops_movers_AlignByResidueSelectorMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace fold_from_loops {
namespace movers {

class AlignByResidueSelectorMover;

typedef utility::pointer::shared_ptr< AlignByResidueSelectorMover > AlignByResidueSelectorMoverOP;
typedef utility::pointer::shared_ptr< AlignByResidueSelectorMover const > AlignByResidueSelectorMoverCOP;

}
} //protocols
} //fold_from_loops

#endif //INCLUDED_protocols_fold_from_loops_AlignByResidueSelectorMover_fwd_hh
