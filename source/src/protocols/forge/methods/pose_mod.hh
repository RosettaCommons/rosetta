// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/methods/pose_mod.hh
/// @brief methods for Pose modifications
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_methods_pose_mod_hh
#define INCLUDED_protocols_forge_methods_pose_mod_hh

// type headers
#include <core/types.hh>

// project headers

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/AtomID.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pose/variant_util.hh>

namespace protocols {
namespace forge {
namespace methods {


/// @brief grow a series of residues to the left of a position
/// @param[in,out] pose Pose the pose to modify
/// @param[in] anchor the anchor position
/// @param[in] begin iterator that points to the first ResidueOP
/// @param[in] end iterator that points just beyond the last ResidueOP
/// @param[in] correct_terminus re-add lower terminus if found, default true
/// @param[in] use_existing_crd No idealization: place residues with existing
///  coordinates, default false.
/// @return the left endpoint of the new growth or zero if nothing to do
/// @remarks Use a reverse iterator to grow right -> left along vector containers.
template< typename ResidueOPIterator >
core::Size
grow_left_r(
	core::pose::Pose & pose,
	core::Size anchor,
	ResidueOPIterator begin,
	ResidueOPIterator end,
	bool const correct_terminus = true,
	bool const use_existing_crd = false
)
{
	using core::Size;

	if ( begin == end ) {
		return 0;
	}

	// query to see if position we're extending from is a true terminus
	bool const had_lower_terminus = pose.residue( anchor ).is_lower_terminus();

	// grow extension
	Size current_pos = anchor; // tracks the anchor as it moves
	for ( ResidueOPIterator i = begin; i != end; ++i, ++current_pos ) {
		pose.conformation().safely_prepend_polymer_residue_before_seqpos( **i, anchor, !use_existing_crd ); // will remove any terminus
	}

	// add terminus if necessary
	if ( correct_terminus && had_lower_terminus && !pose.residue( anchor ).is_lower_terminus() ) {
		core::pose::add_lower_terminus_type_to_pose_residue( pose, anchor );
	}

	return anchor; // the left endpoint of the new growth
}


/// @brief grow a series of residues to the left of a position
/// @param[in,out] pose Pose the pose to modify
/// @param[in] anchor the anchor position
/// @param[in] begin iterator that points to the first ResidueTypeOP
/// @param[in] end iterator that points just beyond the last ResidueTypeOP
/// @param[in] correct_terminus re-add lower terminus if found, default true
/// @return the left endpoint of the new growth or zero if nothing to do
/// @remarks Use a reverse iterator to grow right -> left along vector containers.
template< typename ResidueTypeOPIterator >
core::Size
grow_left_rtype(
	core::pose::Pose & pose,
	core::Size anchor,
	ResidueTypeOPIterator begin,
	ResidueTypeOPIterator end,
	bool const correct_terminus = true
)
{
	using core::conformation::ResidueOPs;
	using core::conformation::ResidueFactory;

	// create residues
	ResidueOPs residues;
	for ( ResidueTypeOPIterator i = begin; i != end; ++i ) {
		residues.push_back( ResidueFactory::create_residue( **i ) );
	}

	return grow_left_r( pose, anchor, residues.begin(), residues.end(), correct_terminus, false );
}


/// @brief grow a series of residues to the right of a position
/// @param[in,out] pose Pose the pose to modify
/// @param[in] anchor the anchor position, can be 0 if operating on empty
///  Pose
/// @param[in] begin iterator that points to the first ResidueOP
/// @param[in] end iterator that points just beyond the last ResidueOP
/// @param[in] correct_terminus re-add upper terminus if found, default true
/// @param[in] use_existing_crd No idealization: place residues with existing
///  coordinates, default false.
/// @return the right endpoint of the new growth or zero if nothing to do
template< typename ResidueOPIterator >
core::Size
grow_right_r(
	core::pose::Pose & pose,
	core::Size anchor,
	ResidueOPIterator begin,
	ResidueOPIterator end,
	bool const correct_terminus = true,
	bool const use_existing_crd = false
)
{
	using core::Size;

	if ( begin == end ) {
		return 0;
	}

	// query to see if position we're extending from is a true terminus
	bool const had_upper_terminus = anchor > 0 ? pose.residue( anchor ).is_upper_terminus() : false;

	// grow extension
	Size current_pos = anchor; // tracks the right endpoint of the new growth

	ResidueOPIterator i = begin;
	if ( anchor == 0 ) { // special case where pose is empty
		pose.append_residue_by_bond( **i );
		++i;
		++current_pos;
	}

	for ( ; i != end; ++i, ++current_pos ) {
		pose.conformation().safely_append_polymer_residue_after_seqpos( **i, current_pos, !use_existing_crd ); // will remove any terminus
	}

	// add terminus if necessary
	if ( correct_terminus && had_upper_terminus && !pose.residue( current_pos ).is_upper_terminus() ) {
		core::pose::add_upper_terminus_type_to_pose_residue( pose, current_pos );
	}

	return current_pos; // the right endpoint of the new growth
}


/// @brief grow a series of residues to the right of a position
/// @param[in,out] pose Pose the pose to modify
/// @param[in] anchor the anchor position, can be 0 if operating on empty
///  Pose
/// @param[in] begin iterator that points to the first ResidueTypeOP
/// @param[in] end iterator that points just beyond the last ResidueTypeOP
/// @param[in] correct_terminus re-add upper terminus if found, default true
/// @return the right endpoint of the new growth or zero if nothing to do
template< typename ResidueTypeOPIterator >
core::Size
grow_right_rtype(
	core::pose::Pose & pose,
	core::Size anchor,
	ResidueTypeOPIterator begin,
	ResidueTypeOPIterator end,
	bool const correct_terminus = true
)
{
	using core::conformation::ResidueOPs;
	using core::conformation::ResidueFactory;

	// create residues
	ResidueOPs residues;
	for ( ResidueTypeOPIterator i = begin; i != end; ++i ) {
		residues.push_back( ResidueFactory::create_residue( **i ) );
	}

	return grow_right_r( pose, anchor, residues.begin(), residues.end(), correct_terminus, false );
}


/// @brief add cutpoint variants at a specific position
/// @param[in,out] Pose to modify
/// @param[in] position at which to add cutpoint variants
/// @return true if cutpoint variants added, false if position not a topological cutpoint
bool
add_cutpoint_variants(
	core::pose::Pose & pose,
	core::Size const pos
);


/// @brief remove cutpoint variants at a specific position
/// @param[in,out] Pose to modify
/// @param[in] position at which to remove cutpoint variants
/// @return true if cutpoint variants removed, false if no cutpoint variants found
///  or position not a topological cutpoint
bool
remove_cutpoint_variants(
	core::pose::Pose & pose,
	core::Size const pos
);


/// @brief restore residues (i.e. sidechains)
/// @param[in] old2new map indicating residues to be transferred and
///  the mapping from archive_pose position -> pose position
/// @param[in] archive_pose the original Pose to take residues from
/// @param[out] pose the altered Pose that needs residue restoration
void
restore_residues(
	std::map< core::Size, core::Size > const & old2new,
	core::pose::Pose & archive_pose,
	core::pose::Pose & pose
);


/// @brief restore residues (i.e. sidechains)
/// @param[in] archive_pose the original Pose to take residues from
/// @param[out] pose the altered Pose that needs residue restoration
/// @remarks length between two poses must be equal
void
restore_residues(
	core::pose::Pose & archive_pose,
	core::pose::Pose & pose
);


} // methods
} // forge
} // protocols


#endif /* INCLUDED_protocols_forge_methods_pose_mod_HH */
