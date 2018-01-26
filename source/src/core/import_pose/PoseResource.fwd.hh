// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/PoseResource.fwd.hh
/// @brief  Declaration for the bitwise-constant class that holds a load-once Pose in memory
///         and protects it so that it can be shared safely between threads.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_import_pose_PoseResource_FWD_HH
#define INCLUDED_core_import_pose_PoseResource_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace import_pose {

class PoseResource;
typedef utility::pointer::shared_ptr< PoseResource > PoseResourceOP;
typedef utility::pointer::shared_ptr< PoseResource const > PoseResourceCOP;

} // namespace import_pose
} // namespace core

#endif //INCLUDED_core_import_pose_PoseResource_FWD_HH

