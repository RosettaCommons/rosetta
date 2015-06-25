// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/pose/reference_pose/ReferencePoseSet.cc
/// @brief  Headers for ReferencePoseSet, a class for holding sets of reference poses.
/// @details Reference poses are a means of storing information about the state of a pose
/// at one point in a protocol and retrieving it later.  The primary usage case is if a
/// pose is going to have an unknown number of residues inserted into it, but certain movers
/// must be set up with reference to residue indices (that might change).  By creating a
/// reference pose, setting up movers with respect to the indices of the reference pose, and
/// tracking how residue indices in the modified pose correspond to residue indices in the
/// reference pose, movers can figure out which residues they actually should be operating on.
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory.

// Unit header
#include <core/pose/reference_pose/ReferencePoseSet.hh>

// Package headers
#include <core/pose/reference_pose/ReferencePose.hh>

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

			static thread_local basic::Tracer TR( "core.pose.reference_pose.ReferencePoseSet" );

			/// @brief Constructor.
			///
			ReferencePoseSet::ReferencePoseSet() :
				reference_pose_map_()
			{
			}

			/// @brief Copy constructor.
			///
			ReferencePoseSet::ReferencePoseSet( ReferencePoseSet const & src ) :
				utility::pointer::ReferenceCount(),
				utility::pointer::enable_shared_from_this< ReferencePoseSet >(),
				reference_pose_map_()
			{
				//Loop over and clone all elements in the map:
				for(std::map<std::string,ReferencePoseOP>::const_iterator it=src.reference_pose_map_.begin(); it!=src.reference_pose_map_.end(); ++it) {
					reference_pose_map_[it->first] = it->second->clone(); //Create a clone of each ReferencePose object.
				}
			}

			/// @brief Destructor.
			///
			ReferencePoseSet::~ReferencePoseSet() {}


			/// @brief Make a copy of this ReferencePoseSet object (allocate actual memory for it)
			/// and return an owning pointer to the copy.
			ReferencePoseSetOP
			ReferencePoseSet::clone() const
			{
				return ReferencePoseSetOP( new ReferencePoseSet( *this ) );
			}
			
			/// @brief Given a reference pose, create a new ReferencePose object and initialize it
			/// based on the reference pose, then add it to the current map.
			/// @details After this operation, the map holds a new ReferencePoseOP pointing to a
			/// very uninteresting ReferencePose object that maps a list of residues onto themselves,
			/// accessible by the key string.
			void ReferencePoseSet::add_and_initialize_reference_pose(
				std::string const &key_string,
				core::pose::Pose const &pose
			) {
				runtime_assert_string_msg(
					key_string != "",
					"Error in core::pose::reference_pose::ReferencePoseSet::add_and_initialize_reference_pose():  The name provided for this reference pose must not be an empty string."
				);
				runtime_assert_string_msg(
					reference_pose_map_.count(key_string)==0,
					"Error in core::pose::reference_pose::ReferencePoseSet::add_and_initialize_reference_pose(): The name provided for this reference pose has already been assigned to another reference pose."
				);
				reference_pose_map_[key_string] = ReferencePoseOP( new ReferencePose() );
				reference_pose_map_[key_string]->initialize_residue_map_from_pose( pose );
				return;
			}
			
			/// @brief Get the corresponding residue in the current pose to a residue in a reference pose.
			/// @details Should return 0 if no corresponding residue exists in the current pose, and throw an
			/// error if the residue index provided as a key did not exist in the reference pose.
			core::Size ReferencePoseSet::corresponding_residue_in_current(
				core::Size const res_index_in_ref,
				std::string const &ref_pose_key_string
			) const {
				runtime_assert_string_msg(
					reference_pose_map_.count(ref_pose_key_string)!=0,
					"Error in core::pose::reference_pose::ReferencePoseSet::corresponding_residue_in_current(): A reference pose with the name " + ref_pose_key_string + " does not exist in the reference pose map."
				);
				return reference_pose_map_.at(ref_pose_key_string)->corresponding_residue_in_current( res_index_in_ref );
			}
			
			/// @brief Get the number of ReferencePose objects in this ReferencePoseSet.
			///					
			core::Size ReferencePoseSet::n_reference_pose() const {
				return reference_pose_map_.size();
			}
			
			/// @brief Returns true if and only if there are no ReferencePose objects in this ReferencePoseSet.
			///					
			bool ReferencePoseSet::empty() const {
				return reference_pose_map_.empty();
			}
			
			/// @brief Find all mappings to indices in the new pose after seqpos in all ReferencePose objects, and increment them by 1.
			/// @details If there is no ReferencePose object, do nothing.
			void ReferencePoseSet::increment_reference_pose_mapping_after_seqpos( core::Size const seqpos ) {
				if(empty()) return; //Do nothing if the map is empty.
				//Iterate over all elements of the ReferencePose map:
				for (std::map< std::string, ReferencePoseOP >::iterator it=reference_pose_map_.begin(); it!=reference_pose_map_.end(); ++it) {
					debug_assert( it->second ); //The reference pose owning pointer should always point to something, but let's double-check that.
					it->second->increment_reference_pose_mapping_after_seqpos( seqpos );
				}
				return;
			}

			/// @brief Find all mappings to indices in the new pose after seqpos in all ReferencePose objects, and decrement them by 1.
			/// @details If there is no ReferencePose object, do nothing.
			void ReferencePoseSet::decrement_reference_pose_mapping_after_seqpos( core::Size const seqpos ) {
				if(empty()) return; //Do nothing if the map is empty.
				//Iterate over all elements of the ReferencePose map:
				for (std::map< std::string, ReferencePoseOP >::iterator it=reference_pose_map_.begin(); it!=reference_pose_map_.end(); ++it) {
					debug_assert( it->second ); //The reference pose owning pointer should always point to something, but let's double-check that.
					it->second->decrement_reference_pose_mapping_after_seqpos( seqpos );
				}
				return;
			} 

			/// @brief Find all mappings to indices in the new pose to seqpos in all ReferencePose objects, and set them to point to residue 0 (deletion signal).
			/// @details If there is no ReferencePose object, do nothing.
			void ReferencePoseSet::zero_reference_pose_mapping_at_seqpos( core::Size const seqpos ) {
				if(empty()) return; //Do nothing if the map is empty.
				//Iterate over all elements of the ReferencePose map:
				for (std::map< std::string, ReferencePoseOP >::iterator it=reference_pose_map_.begin(); it!=reference_pose_map_.end(); ++it) {
					debug_assert( it->second ); //The reference pose owning pointer should always point to something, but let's double-check that.
					it->second->zero_reference_pose_mapping_at_seqpos( seqpos );
				}
				return;
			}

		} // namespace reference_pose
	} // namespace pose
} // namespace core

