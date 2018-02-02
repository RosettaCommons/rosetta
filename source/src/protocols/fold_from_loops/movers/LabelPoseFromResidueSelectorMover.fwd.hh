// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/LabelPoseFromResidueSelectorMover.fwd.hh
/// @brief Adds labels to residues selected through a ResidueSelector
/// @author Jaume Bonet (jaume.bonet@gmail.com)


#ifndef INCLUDED_protocols_fold_from_loops_movers_LabelPoseFromResidueSelectorMover_fwd_hh
#define INCLUDED_protocols_fold_from_loops_movers_LabelPoseFromResidueSelectorMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace fold_from_loops {
namespace movers {

class LabelPoseFromResidueSelectorMover;

typedef utility::pointer::shared_ptr< LabelPoseFromResidueSelectorMover > LabelPoseFromResidueSelectorMoverOP;
typedef utility::pointer::shared_ptr< LabelPoseFromResidueSelectorMover const > LabelPoseFromResidueSelectorMoverCOP;

}
} //protocols
} //fold_from_loops

#endif //INCLUDED_protocols_fold_from_loops_LabelPoseFromResidueSelectorMover_fwd_hh
