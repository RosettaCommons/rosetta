// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/Pose.fwd.hh
/// @brief  Pose forward declarations header
/// @author Phil Bradley

#ifndef INCLUDED_core_pose_Pose_fwd_hh
#define INCLUDED_core_pose_Pose_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace core {
namespace pose {

// Forward
class Pose;

using PoseOP  = utility::pointer::shared_ptr< Pose > ;
using PoseCOP = utility::pointer::shared_ptr< Pose const >;

using PoseAP  = utility::pointer::weak_ptr< Pose >;
using PoseCAP = utility::pointer::weak_ptr< Pose const >;

using PoseUP  = std::unique_ptr< Pose >;
using PoseCUP = std::unique_ptr< Pose const >;

using PoseOPs  = utility::vector1< PoseOP >;
using PoseCOPs = utility::vector1< PoseCOP >;

} // namespace pose
} // namespace core

#endif // INCLUDED_core_pose_Pose_FWD_HH
