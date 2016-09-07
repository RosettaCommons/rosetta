// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/PoseFromPDBLoader.hh
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_core_import_pose_PoseFromPDBLoader_HH
#define INCLUDED_core_import_pose_PoseFromPDBLoader_HH

//unit headers
#include <core/import_pose/PoseFromPDBLoader.fwd.hh>

//project headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>


//utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <istream>

namespace core {
namespace import_pose {

class PoseFromPDBLoader : public basic::resource_manager::ResourceLoader
{
public:
	PoseFromPDBLoader();
	~PoseFromPDBLoader() override;

	/// @brief Returns an owning pointer to a import_pose_options object
	/// which is constructed from the given input stream (istream) which in tern
	/// originates from a particular data source (given by the name input_tag)

	utility::pointer::ReferenceCountOP
	create_resource(
		basic::resource_manager::ResourceOptions const & options,
		basic::resource_manager::LocatorID const & locator_id,
		std::istream & istream
	) const override;


	basic::resource_manager::ResourceOptionsOP
	default_options() const override;

};

} // namespace import_pose
} // namespace core

#endif //INCLUDED_core_import_pose_PoseFromPDBLoader_HH

