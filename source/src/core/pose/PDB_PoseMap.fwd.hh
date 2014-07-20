// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/PDB_PoseMap.fwd.hh
/// @brief  fwd declaration of classes defined in pose/PDB_PoseMap
/// @author Steven Lewis


#ifndef INCLUDED_core_pose_PDB_PoseMap_fwd_hh
#define INCLUDED_core_pose_PDB_PoseMap_fwd_hh


// Unit headers

// Package headers

// Project headers

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ Headers
namespace core {
namespace pose {

class PDB_PoseMap;

typedef utility::pointer::owning_ptr< PDB_PoseMap > PDB_PoseMapOP;
typedef utility::pointer::owning_ptr< PDB_PoseMap const > PDB_PoseMapCOP;

} // namespace pose
} // namespace core

#ifdef USEBOOSTSERIALIZE
#include <boost/serialization/map.hpp>
#endif


#endif
