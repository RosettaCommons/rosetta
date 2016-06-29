// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/movers/RotateSegmentMover.fwd.hh
/// @brief Rotates a segment in the pose
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_movers_RotateSegmentMover_fwd_hh
#define INCLUDED_protocols_denovo_design_movers_RotateSegmentMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace denovo_design {
namespace movers {

class RotateSegmentMover;

typedef utility::pointer::shared_ptr< RotateSegmentMover > RotateSegmentMoverOP;
typedef utility::pointer::shared_ptr< RotateSegmentMover const > RotateSegmentMoverCOP;

} //protocols
} //denovo_design
} //movers

#endif //INCLUDED_protocols_denovo_design_movers_RotateSegmentMover_fwd_hh

