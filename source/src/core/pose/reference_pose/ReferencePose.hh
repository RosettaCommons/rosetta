// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/reference_pose/ReferencePose.hh
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


#ifndef INCLUDED_core_pose_reference_pose_ReferencePose_hh
#define INCLUDED_core_pose_reference_pose_ReferencePose_hh


// Unit headers
#include <core/pose/reference_pose/ReferencePose.fwd.hh>

// Package headers
#include <core/pose/Pose.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers

// C++ headers
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace reference_pose {

/// @brief  ReferencePose class, used to store sets of ReferencePose data for tracking how a pose changes over the course of a protocol.
///
class ReferencePose : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< ReferencePose >
{
public:

	/// @brief Constructor
	///
	ReferencePose();

	/// @brief Copy constructor
	///
	ReferencePose( ReferencePose const & src );

	/// @brief Destructor
	/// @details Note that Rosetta destructors are generally not called.
	~ReferencePose();

	/// @brief Make a copy of this ReferencePose object (allocate actual memory for it)
	/// and return an owning pointer to the copy.
	virtual
	ReferencePoseOP clone() const;

	/// @brief Get const owning pointer (e.g. from a reference).
	///
	inline ReferencePoseCOP get_self_ptr() const { return shared_from_this(); }

	/// @brief Get owning pointer (e.g. from a reference).
	///
	inline ReferencePoseOP get_self_ptr() { return shared_from_this(); }

	/// @brief Get const access (weak) pointer (e.g. from a reference).
	///
	inline ReferencePoseCAP get_self_weak_ptr() const { return ReferencePoseCAP( shared_from_this() ); }

	/// @brief Get access (weak) pointer (e.g. from a reference).
	///
	inline ReferencePoseAP get_self_weak_ptr() { return ReferencePoseAP( shared_from_this() ); }

	/********************************************************************************
	GETTERS
	*********************************************************************************/

	/// @brief Returns true if this ReferencePose object stores a map of
	/// (old residue indices->new residue indices), false otherwise.
	bool stores_residue_map() const { return stores_residue_map_; }

	/********************************************************************************
	SETTERS
	*********************************************************************************/

	/// @brief Initializes the residue_map_ based on an input reference pose.
	/// @details After this operation, the residue_map_ is a very uninteresting map
	/// that maps every residue index in the pose onto itself.  This also sets the
	/// stores_residue_map_ bit to true.
	void initialize_residue_map_from_pose( core::pose::Pose const &pose );

	/// @brief Find all mappings to indices in the new pose after seqpos, and increment them by 1.
	/// @details If there is no ReferencePose object, do nothing.
	void increment_reference_pose_mapping_after_seqpos( core::Size const seqpos );

	/// @brief Find all mappings to indices in the new pose after seqpos, and decrement them by 1.
	/// @details If there is no ReferencePose object, do nothing.
	void decrement_reference_pose_mapping_after_seqpos( core::Size const seqpos );

	/// @brief Find all mappings to indices in the new pose to seqpos, and set them to point to residue 0 (deletion signal).
	/// @details If there is no ReferencePose object, do nothing.
	void zero_reference_pose_mapping_at_seqpos( core::Size const seqpos );

	/********************************************************************************
	GETTERS
	*********************************************************************************/

	/// @brief Get the residue in the current pose corresponding to a particular residue in the
	/// reference pose.
	/// @details Should return 0 if there is no corresponding residue in the current pose (e.g. if the
	/// residue was deleted.)  Throws an error if there was no residue with the given index in the
	/// reference pose.
	core::Size corresponding_residue_in_current( core::Size const res_in_ref ) const {
		runtime_assert_string_msg( residue_map_.count(res_in_ref)!=0, "Error in core::pose::reference_pose::ReferencePose::corresponding_residue_in_current(): The given residue index did not exist in the reference pose." );
		return residue_map_.at(res_in_ref);
	}

private:

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

	/// @brief Does this ReferencePose store a residue map of (old residue indices -> new residue indices)?
	/// @details Default false.
	bool stores_residue_map_;

	/// @brief Mapping of reference pose residue indices (key) onto new pose
	/// residue indices (mapped values).
	std::map <core::Size, core::Size> residue_map_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class ReferencePose

} // namespace reference_pose
} // namespace pose
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_reference_pose_ReferencePose )
#endif // SERIALIZATION


#endif
