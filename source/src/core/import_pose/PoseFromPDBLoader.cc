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

#include <basic/resource_manager/ResourceManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/PoseResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>

//utility headers
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

// numeric headers

static basic::Tracer TR( "core.import_pose.PoseFromPDBLoader" );

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
		throw CREATE_EXCEPTION(utility::excn::Exception,  "PoseFromPDBLoader expected to be given a ImportPoseOptions object, " \
			"but was given a non-ImportPoseOptions object of type '" + options.type() + "', which has the name '" + options.name() + "'." );
	}
	pose::PoseOP pose( new pose::Pose() );
	//////////////////////////////////////////////////////////////////////////
	// This inappropriately specific code is needed to cross load ResidueTypes
	// into the pose loader
	//
	basic::resource_manager::ResourceManager * resource_manager( basic::resource_manager::ResourceManager::get_instance());
	if ( resource_manager != nullptr && resource_manager->has_resource_with_description( "residue" ) ) {
		basic::resource_manager::ResourceOP residue_resource( resource_manager->get_resource( "residue" ) );
		core::chemical::ResidueTypeOP new_residue(utility::pointer::dynamic_pointer_cast< core::chemical::ResidueType > ( residue_resource ));
		debug_assert( new_residue );

		core::chemical::ResidueTypeSet const & orig_rts( *pose->residue_type_set_for_pose( new_residue->mode() ) );
		if ( ! orig_rts.has_name( new_residue->name() ) ) {
			TR << "loading residue " << new_residue->name() << " into " << new_residue->mode() <<" mode residue_type_set" <<std::endl;
			core::chemical::PoseResidueTypeSetOP new_rts( pose->conformation().modifiable_residue_type_set_for_conf( new_residue->mode() ) );
			// TODO: Should this be patchable?
			new_rts->add_unpatchable_residue_type(new_residue);
			pose->conformation().reset_residue_type_set_for_conf( new_rts );
		}
	}
	///////////////////////////////////////////////////////////////////////////
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
