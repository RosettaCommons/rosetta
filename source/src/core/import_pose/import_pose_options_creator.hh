// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/import_pose_options_creator.hh
/// @brief
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_core_import_pose_import_pose_options_creator_HH
#define INCLUDED_core_import_pose_import_pose_options_creator_HH

//unit headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>

namespace core {
namespace import_pose {

class ImportPoseOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
{
public:
	ImportPoseOptionsCreator();
	virtual ~ImportPoseOptionsCreator();

	virtual basic::resource_manager::ResourceOptionsOP create_options() const;
	virtual std::string options_type() const;

};

} // namespace import_pose
} // namespace core

#endif //INCLUDED_core_import_pose_import_pose_options_creator_HH
