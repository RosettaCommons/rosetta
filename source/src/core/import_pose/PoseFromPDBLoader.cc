// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/PoseFromPDBLoader.cc
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

//unit headers
#include <core/import_pose/PoseFromPDBLoader.hh>
#include <core/import_pose/PoseFromPDBLoaderCreator.hh>

//package headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>

//utility headers
#include <utility/excn/Exceptions.hh>

// numeric headers

namespace core {
namespace import_pose {

PoseFromPDBLoader::PoseFromPDBLoader() {}
PoseFromPDBLoader::~PoseFromPDBLoader() = default;

utility::pointer::ReferenceCountOP
PoseFromPDBLoader::create_resource(
	basic::resource_manager::ResourceOptions const & options,
	basic::resource_manager::LocatorID const & locator_id,
	std::istream & istream
) const
{
	ImportPoseOptions const * pose_opts_ptr = dynamic_cast< ImportPoseOptions const * > ( &options );
	if ( ! pose_opts_ptr ) {
		throw utility::excn::EXCN_Msg_Exception( "PoseFromPDBLoader expected to be given a ImportPoseOptions object, " \
			"but was given a non-ImportPoseOptions object of type '" + options.type() + "', which has the name '" + options.name() + "'." );
	}
	pose::PoseOP pose( new pose::Pose() );
	pose_from_pdb_stream( *pose, istream, locator_id, *pose_opts_ptr );
	return pose;
}

basic::resource_manager::ResourceOptionsOP
PoseFromPDBLoader::default_options() const
{
	return basic::resource_manager::ResourceOptionsOP( new ImportPoseOptions() );
}

basic::resource_manager::ResourceLoaderOP PoseFromPDBLoaderCreator::create_resource_loader() const
{
	return basic::resource_manager::ResourceLoaderOP( new PoseFromPDBLoader() );
}

std::string PoseFromPDBLoaderCreator::loader_type() const
{
	return "PoseFromPDB";
}

} // namespace import_pose
} // namespace core
