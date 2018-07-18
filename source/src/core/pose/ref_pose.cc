// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/util.cc
/// @brief  Pose class utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Rhiju Das, Steven Lewis, Vikram K. Mulligan


// Unit header
#include <core/pose/ref_pose.hh>

// Package headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/Tag.hh>

// C/C++ headers
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <boost/functional/hash.hpp>



namespace core {
namespace pose {

static basic::Tracer TR( "core.pose.ref_pose" );


void
attributes_for_saved_reference_pose(
	utility::tag::AttributeList & attributes,
	std::string const & attribute_name)
{
	attributes_for_saved_reference_pose_w_description( attributes, "", attribute_name );
}

void
attributes_for_saved_reference_pose_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description,
	std::string const & attribute_name)
{
	attributes + utility::tag::XMLSchemaAttribute( attribute_name, utility::tag::xs_string,
		( description == "" ? "Name of reference pose to use" : description ) );
}

core::pose::PoseOP
saved_reference_pose( utility::tag::TagCOP const in_tag, basic::datacache::DataMap & data_map, std::string const & tag_name ){

	if ( in_tag->hasOption(tag_name) ) {
		core::pose::PoseOP refpose(nullptr);
		std::string refpose_name(in_tag->getOption<std::string>( tag_name) );
		TR<<"Loading PDB: "<<refpose_name<<std::endl;

		if ( !data_map.has("spm_ref_poses",refpose_name) ) {
			refpose = core::pose::PoseOP( new core::pose::Pose() );
			data_map.add("spm_ref_poses",refpose_name,refpose );
		} else refpose = data_map.get_ptr<core::pose::Pose>("spm_ref_poses",refpose_name );

		return refpose;
	} else std::cerr << "WARNING: saved_reference_pose function called even though tag has no " + tag_name + " entry. something's unclean somewhere." << std::endl;
	return nullptr;
}

core::pose::PoseCOP
saved_native_pose( basic::datacache::DataMap & data_map, std::string const & resource_str ){

	if ( data_map.has_resource( resource_str) ) {
		core::pose::PoseCOP native = data_map.get_resource<core::pose::Pose const>( resource_str );
		return native;
	} else {
		TR << "Native Pose Resource not found in datamap: " << resource_str << std::endl;
		return nullptr;
	}
}

} // pose
} // core
