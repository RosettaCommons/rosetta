// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/PoseFromPoseResourceMover.fwd.hh
/// @brief Mover that shuttles a Pose from a PoseResource into the DataMap when its parse_my_tag method is invoked
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_PoseFromPoseResourceMover_fwd_hh
#define INCLUDED_protocols_simple_moves_PoseFromPoseResourceMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace simple_moves {

class PoseFromPoseResourceMover;

typedef utility::pointer::shared_ptr< PoseFromPoseResourceMover > PoseFromPoseResourceMoverOP;
typedef utility::pointer::shared_ptr< PoseFromPoseResourceMover const > PoseFromPoseResourceMoverCOP;

} //protocols
} //simple_moves

#endif //INCLUDED_protocols_simple_moves_PoseFromPoseResourceMover_fwd_hh
