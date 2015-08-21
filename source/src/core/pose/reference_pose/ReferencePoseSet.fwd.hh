// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/reference_pose/ReferencePoseSet.fwd.hh
/// @brief  Forward declarations for ReferencePoseSet, a class for holding sets of reference poses.
/// @details Reference poses are a means of storing information about the state of a pose
/// at one point in a protocol and retrieving it later.  The primary usage case is if a
/// pose is going to have an unknown number of residues inserted into it, but certain movers
/// must be set up with reference to residue indices (that might change).  By creating a
/// reference pose, setting up movers with respect to the indices of the reference pose, and
/// tracking how residue indices in the modified pose correspond to residue indices in the
/// reference pose, movers can figure out which residues they actually should be operating on.
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory.


#ifndef INCLUDED_core_pose_reference_pose_ReferencePoseSet_fwd_hh
#define INCLUDED_core_pose_reference_pose_ReferencePoseSet_fwd_hh

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace core {
namespace pose {
namespace reference_pose {

class ReferencePoseSet;

typedef  utility::pointer::weak_ptr< ReferencePoseSet >  ReferencePoseSetAP;
typedef  utility::pointer::weak_ptr< ReferencePoseSet const >  ReferencePoseSetCAP;
typedef  utility::pointer::shared_ptr< ReferencePoseSet >  ReferencePoseSetOP;
typedef  utility::pointer::shared_ptr< ReferencePoseSet const >  ReferencePoseSetCOP;

typedef  utility::vector1< ReferencePoseSetOP >  ReferencePoseSetOPs;
typedef  utility::vector1< ReferencePoseSetCOP >  ReferencePoseSetCOPs;
typedef  utility::vector1< ReferencePoseSetCAP >  ReferencePoseSetCAPs;

} // namespace reference_pose
} // namespace pose
} // namespace core

#endif // INCLUDED_core_pose_reference_pose_ReferencePoseSet_fwd_hh
