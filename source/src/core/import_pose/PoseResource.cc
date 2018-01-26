// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/PoseResource.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

//unit headers
#include <core/import_pose/PoseResource.hh>

//package headers
#include <core/pose/Pose.hh>

namespace core {
namespace import_pose {

PoseResource::PoseResource( core::pose::PoseCOP pose ) : pose_( pose ) {}
PoseResource::~PoseResource() = default;

void PoseResource::pose( core::pose::PoseCOP pose_to_hold )
{
	pose_ = pose_to_hold;
}

core::pose::PoseOP
PoseResource::pose_deep_copy() const
{
	core::pose::PoseOP copy( new core::pose::Pose );
	copy->detached_copy( *pose_ );
	return copy;
}


} // namespace import_pose
} // namespace core
