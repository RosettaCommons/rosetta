// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/PDBPoseMap.fwd.hh
/// @brief  fwd declaration of classes defined in pose/PDBPoseMap
/// @author Steven Lewis


#ifndef INCLUDED_core_pose_PDBPoseMap_fwd_hh
#define INCLUDED_core_pose_PDBPoseMap_fwd_hh


// Unit headers

// Package headers

// Project headers

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ Headers
namespace core {
namespace pose {

class PDBPoseMap;

typedef utility::pointer::shared_ptr< PDBPoseMap > PDBPoseMapOP;
typedef utility::pointer::shared_ptr< PDBPoseMap const > PDBPoseMapCOP;

} // namespace pose
} // namespace core


#endif
