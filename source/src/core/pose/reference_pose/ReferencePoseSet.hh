// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/reference_pose/ReferencePoseSet.hh
/// @brief  Headers for ReferencePoseSet, a class for holding sets of reference poses.
/// @details Reference poses are a means of storing information about the state of a pose
/// at one point in a protocol and retrieving it later.  The primary usage case is if a
/// pose is going to have an unknown number of residues inserted into it, but certain movers
/// must be set up with reference to residue indices (that might change).  By creating a
/// reference pose, setting up movers with respect to the indices of the reference pose, and
/// tracking how residue indices in the modified pose correspond to residue indices in the
/// reference pose, movers can figure out which residues they actually should be operating on.
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory.


#ifndef INCLUDED_core_pose_reference_pose_ReferencePoseSet_hh
#define INCLUDED_core_pose_reference_pose_ReferencePoseSet_hh


// Unit headers
#include <core/pose/reference_pose/ReferencePoseSet.fwd.hh>

// Package headers
#include <core/pose/reference_pose/ReferencePose.fwd.hh>
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

/// @brief  ReferencePoseSet class, used to store sets of ReferencePose data for tracking how a pose changes over the course of a protocol.
///
class ReferencePoseSet : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< ReferencePoseSet >
{
public:

	/// @brief Constructor
	///
	ReferencePoseSet();

	/// @brief Copy constructor
	///
	ReferencePoseSet( ReferencePoseSet const & src );

	/// @brief Destructor
	/// @details Note that Rosetta destructors are generally not called.
	~ReferencePoseSet();

	/// @brief Make a copy of this ReferencePoseSet object (allocate actual memory for it)
	/// and return an owning pointer to the copy.
	virtual
	ReferencePoseSetOP clone() const;

	/// @brief Get const owning pointer (e.g. from a reference).
	///
	inline ReferencePoseSetCOP get_self_ptr() const { return shared_from_this(); }

	/// @brief Get owning pointer (e.g. from a reference).
	///
	inline ReferencePoseSetOP get_self_ptr() { return shared_from_this(); }

	/// @brief Get const access (weak) pointer (e.g. from a reference).
	///
	inline ReferencePoseSetCAP get_self_weak_ptr() const { return ReferencePoseSetCAP( shared_from_this() ); }

	/// @brief Get access (weak) pointer (e.g. from a reference).
	///
	inline ReferencePoseSetAP get_self_weak_ptr() { return ReferencePoseSetAP( shared_from_this() ); }

	/// @brief Given a reference pose, create a new ReferencePose object and initialize it
	/// based on the reference pose, then add it to the current map.
	/// @details After this operation, the map holds a new ReferencePoseOP pointing to a
	/// very uninteresting ReferencePose object that maps a list of residues onto themselves,
	/// accessible by the key string.
	/// If override_current is set to true, will overrwrite any reference pose with the key_string.

	void add_and_initialize_reference_pose( std::string const &key_string, core::pose::Pose const &pose, bool override_current = false );

	/// @brief Get an owning pointer to a reference pose stored in the reference pose map, given a key string.
	/// @details An error is thrown if the key string doesn't associate with anything (i.e. if the reference pose
	/// doesn't exist).
	ReferencePoseOP reference_pose( std::string const &key_string ) {
		runtime_assert_string_msg( reference_pose_map_.count( key_string )!=0,
			"Error in core::pose::reference_pose::ReferencePoseSet::reference_pose(): A reference pose with the name " + key_string + " does not exist in the reference pose map." );
		return reference_pose_map_[key_string];
	}

	/// @brief Get a const owning pointer to a reference pose stored in the reference pose map, given a key string.
	/// @details An error is thrown if the key string doesn't associate with anything (i.e. if the reference pose
	/// doesn't exist).
	ReferencePoseCOP reference_pose_cop( std::string const &key_string ) const {
		runtime_assert_string_msg( reference_pose_map_.count( key_string )!=0,
			"Error in core::pose::reference_pose::ReferencePoseSet::reference_pose(): A reference pose with the name " + key_string + " does not exist in the reference pose map." );
		return ReferencePoseCOP( reference_pose_map_.at( key_string ) );
	}

	/// @brief Get the corresponding residue in the current pose to a residue in a reference pose.
	/// @details Should return 0 if no corresponding residue exists in the current pose, and throw an
	/// error if the residue index provided as a key did not exist in the reference pose.
	core::Size corresponding_residue_in_current(
		core::Size const res_index_in_ref,
		std::string const &ref_pose_key_string
	) const;

	/// @brief Get the number of ReferencePose objects in this ReferencePoseSet.
	///
	core::Size n_reference_pose() const;

	/// @brief Returns true if and only if there are no ReferencePose objects in this ReferencePoseSet.
	///
	bool empty() const;

	/// @brief Find all mappings to indices in the new pose after seqpos in all ReferencePose objects, and increment them by 1.
	/// @details If there is no ReferencePose object, do nothing.
	void increment_reference_pose_mapping_after_seqpos( core::Size const seqpos );

	/// @brief Find all mappings to indices in the new pose after seqpos in all ReferencePose objects, and decrement them by 1.
	/// @details If there is no ReferencePose object, do nothing.
	void decrement_reference_pose_mapping_after_seqpos( core::Size const seqpos );

	/// @brief Find all mappings to indices in the new pose to seqpos in all ReferencePose objects, and set them to point to residue 0 (deletion signal).
	/// @details If there is no ReferencePose object, do nothing.
	void zero_reference_pose_mapping_at_seqpos( core::Size const seqpos );


private:

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

	/// @brief Owning pointers to the reference poses in this object, indexed by unique
	/// name strings.
	std::map < std::string, ReferencePoseOP > reference_pose_map_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class ReferencePoseSet

} // namespace reference_pose
} // namespace pose
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_reference_pose_ReferencePoseSet )
#endif // SERIALIZATION


#endif
