// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/reference_pose/ReferencePose.fwd.hh
/// @brief  Forward declarations for ReferencePose, a class for holding information relating
/// the current pose to a reference pose.
/// @details Reference poses are a means of storing information about the state of a pose
/// at one point in a protocol and retrieving it later.  The primary usage case is if a
/// pose is going to have an unknown number of residues inserted into it, but certain movers
/// must be set up with reference to residue indices (that might change).  By creating a
/// reference pose, setting up movers with respect to the indices of the reference pose, and
/// tracking how residue indices in the modified pose correspond to residue indices in the
/// reference pose, movers can figure out which residues they actually should be operating on.
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory.


#ifndef INCLUDED_core_pose_reference_pose_ReferencePose_fwd_hh
#define INCLUDED_core_pose_reference_pose_ReferencePose_fwd_hh

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace core {
namespace pose {
namespace reference_pose {

class ReferencePose;

typedef  utility::pointer::weak_ptr< ReferencePose >  ReferencePoseAP;
typedef  utility::pointer::weak_ptr< ReferencePose const >  ReferencePoseCAP;
typedef  utility::pointer::shared_ptr< ReferencePose >  ReferencePoseOP;
typedef  utility::pointer::shared_ptr< ReferencePose const >  ReferencePoseCOP;

typedef  utility::vector1< ReferencePoseOP >  ReferencePoseOPs;
typedef  utility::vector1< ReferencePoseCOP >  ReferencePoseCOPs;
typedef  utility::vector1< ReferencePoseCAP >  ReferencePoseCAPs;

} // namespace reference_pose
} // namespace pose
} // namespace core

#endif // INCLUDED_core_pose_reference_pose_ReferencePose_fwd_hh
