// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/PoseResourceLoader.hh
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_core_import_pose_PoseResourceLoader_HH
#define INCLUDED_core_import_pose_PoseResourceLoader_HH

//unit headers
#include <core/import_pose/PoseResourceLoader.fwd.hh>

//project headers
#include <basic/resource_manager/ResourceLoader.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

//C++ headers
#include <istream>

namespace core {
namespace import_pose {

class PoseResourceLoader : public basic::resource_manager::ResourceLoader
{
public:
	PoseResourceLoader();
	~PoseResourceLoader() override;

	/// @brief Returns an owning pointer to a import_pose_options object
	/// which is constructed from the given input stream (istream) which in tern
	/// originates from a particular data source (given by the name input_tag)

	basic::resource_manager::ResourceCOP
	create_resource(
		basic::resource_manager::ResourceManager & resource_manager,
		utility::tag::TagCOP resource_tag,
		std::string const & input_id,
		std::istream & input_stream
	) const override;

	static
	std::string
	classname();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


};

} // namespace import_pose
} // namespace core

#endif //INCLUDED_core_import_pose_PoseResourceLoader_HH

