// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/pose/reference_pose/ReferencePose.cc
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

// Unit header
#include <core/pose/reference_pose/ReferencePose.hh>

// Package headers
#include <core/pose/Pose.hh>

// Project headers

// Basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Numeric headers

// Utility Headers
#include <utility/assert.hh>
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <iostream>


namespace core {
namespace pose {
namespace reference_pose {

static thread_local basic::Tracer TR( "core.pose.reference_pose.ReferencePose" );

/// @brief Constructor.
///
ReferencePose::ReferencePose() :
	stores_residue_map_( false ),
	residue_map_()
{
}

/// @brief Copy constructor.
///
ReferencePose::ReferencePose( ReferencePose const & src ) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< ReferencePose >(),
	stores_residue_map_( src.stores_residue_map_ ),
	residue_map_( src.residue_map_ )
{
}

/// @brief Destructor.
///
ReferencePose::~ReferencePose() {}


/// @brief Make a copy of this ReferencePose object (allocate actual memory for it)
/// and return an owning pointer to the copy.
ReferencePoseOP
ReferencePose::clone() const
{
	return ReferencePoseOP( new ReferencePose( *this ) );
}

/// @brief Initializes the residue_map_ based on an input reference pose.
/// @details After this operation, the residue_map_ is a very uninteresting map
/// that maps every residue index in the pose onto itself.  This also sets the
/// stores_residue_map_ bit to true.
void ReferencePose::initialize_residue_map_from_pose( core::pose::Pose const &pose ) {
	residue_map_.clear();
	stores_residue_map_=true;

	if ( pose.n_residue() > 0 ) {
		for ( core::Size ir=1, irmax=pose.n_residue(); ir<=irmax; ++ir ) {
			residue_map_[ir] = ir;
		}
	}
	return;
}

/// @brief Find all mappings to indices in the new pose after seqpos, and increment them by 1.
/// @details If there is no ReferencePose object, do nothing.
void ReferencePose::increment_reference_pose_mapping_after_seqpos( core::Size const seqpos ) {
	if ( residue_map_.empty() ) return; //Do nothing if the residue map is empty.
	//Loop through all elements of the residue_map:
	for ( std::map< core::Size, core::Size >::iterator it=residue_map_.begin(); it!=residue_map_.end(); ++it ) {
		if ( it->second > seqpos ) it->second = (it->second + 1); //If the mapped value is greater than the seqpos, increment it.
	}
	return;
}

/// @brief Find all mappings to indices in the new pose after seqpos, and decrement them by 1.
/// @details If there is no ReferencePose object, do nothing.
void ReferencePose::decrement_reference_pose_mapping_after_seqpos( core::Size const seqpos ) {
	if ( residue_map_.empty() ) return; //Do nothing if the residue map is empty.
	//Loop through all elements of the residue_map:
	for ( std::map< core::Size, core::Size >::iterator it=residue_map_.begin(); it!=residue_map_.end(); ++it ) {
		if ( it->second > seqpos ) it->second = (it->second - 1); //If the mapped value is greater than the seqpos, decrement it.
	}
	return;
}

/// @brief Find all mappings to indices in the new pose to seqpos, and set them to point to residue 0 (deletion signal).
/// @details If there is no ReferencePose object, do nothing.
void ReferencePose::zero_reference_pose_mapping_at_seqpos( core::Size const seqpos ) {
	if ( residue_map_.empty() ) return; //Do nothing if the residue map is empty.
	//Loop through all elements of the residue_map:
	for ( std::map< core::Size, core::Size >::iterator it=residue_map_.begin(); it!=residue_map_.end(); ++it ) {
		if ( it->second == seqpos ) it->second = 0; //If the mapped value matches the seqpos, zero it.
	}
	return;
}

} // namespace reference_pose
} // namespace pose
} // namespace core

