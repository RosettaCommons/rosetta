// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/ref_pose.hh
/// @brief  Reference Pose utilities
/// @author Jared Adolf-Bryfogle
/// @author Moved from protocols to allow ref pose use outside of protocols

#ifndef INCLUDED_core_pose_ref_pose_hh
#define INCLUDED_core_pose_ref_pose_hh

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.tmpl.hh>

// Project headers
#include <core/types.hh>


// Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace pose {


core::pose::PoseOP
saved_reference_pose(
	utility::tag::TagCOP in_tag,
	basic::datacache::DataMap & data_map,
	std::string const & tag_str="reference_name" );

///@brief
/// Retrieve the native pose from the DataMap which is added through the OptionsCollection
/// during RosettaScript parsing as a resource.
///
///@details
/// If the resource is not present, will return a nullptr.
core::pose::PoseCOP
saved_native_pose(
	basic::datacache::DataMap & data_map,
	std::string const & resource_str = "native_pose" );

/// @brief convenience function to access pointers to poses that will be stored
/// in the data map at an arbitrary point during an RS protocol
/// Will look for tag in in_tag variable
void
attributes_for_saved_reference_pose( utility::tag::AttributeList & attributes, std::string const & attribute_name="reference_name" );

/// @brief convenience function to access pointers to poses that will be stored
/// in the data map at an arbitrary point during an RS protocol
/// Will look for tag in in_tag variable
void
attributes_for_saved_reference_pose_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description,
	std::string const & attribute_name="reference_name" );


} // pose
} // core

#endif // INCLUDED_core_pose_ref_pose_HH
