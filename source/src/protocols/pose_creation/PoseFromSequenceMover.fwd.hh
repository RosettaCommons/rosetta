// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_creation/PoseFromSequenceMover.fwd.hh
/// @brief A class for generating a pose from a sequence/fasta
/// @author Dan Farrell (danpf@uw.edu)

#ifndef INCLUDED_protocols_pose_creation_PoseFromSequenceMover_fwd_hh
#define INCLUDED_protocols_pose_creation_PoseFromSequenceMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace pose_creation {

class PoseFromSequenceMover;

using PoseFromSequenceMoverOP = utility::pointer::shared_ptr< PoseFromSequenceMover >;
using PoseFromSequenceMoverCOP = utility::pointer::shared_ptr< PoseFromSequenceMover const >;

} //pose_creation
} //protocols

#endif //INCLUDED_protocols_pose_creation_PoseFromSequenceMover_fwd_hh
