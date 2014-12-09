// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/MiniPose.fwd.hh
/// @brief  MiniPose forward declarations header
/// @author Rhiju Das

#ifndef INCLUDED_core_pose_MiniPose_FWD_HH
#define INCLUDED_core_pose_MiniPose_FWD_HH

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/types.hh>

#include <map>

namespace core {
namespace pose {

// Mapping from these little poses to big pose.
typedef std::map< core::Size, core::Size > ResMap;

// Forward
class MiniPose;

typedef utility::pointer::shared_ptr< MiniPose > MiniPoseOP;
typedef utility::pointer::shared_ptr< MiniPose const > MiniPoseCOP;

typedef utility::pointer::weak_ptr< MiniPose > MiniPoseAP;
typedef utility::pointer::weak_ptr< MiniPose const > MiniPoseCAP;

} // namespace pose
} // namespace core


#endif // INCLUDED_core_pose_MiniPose_FWD_HH
