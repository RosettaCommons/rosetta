// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/surface_docking/SurfaceVectorLoader.hh
/// @brief
/// @author Michael Pacella (mpacella88@gmail.com)

#ifndef INCLUDED_protocols_surface_docking_SurfaceVectorLoader_HH
#define INCLUDED_protocols_surface_docking_SurfaceVectorLoader_HH

//unit headers
#include <protocols/surface_docking/SurfaceVectorLoader.fwd.hh>

// Package headers
#include <protocols/surface_docking/SurfaceParameters.fwd.hh>

//project headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

//C++ headers
#include <istream>

namespace protocols {
namespace surface_docking {

/// @brief %SurfaceVectorLoader constructs a SurfaceParameters instance from data provided by the resource manager
/// @details The SurfaceVectorLoader is given an istream of three points in cartesian space used to configure SurfaceParameters
class SurfaceVectorLoader : public basic::resource_manager::ResourceLoader
{
public:
	SurfaceVectorLoader();
	~SurfaceVectorLoader() override;

	/// @brief Returns a SurfaceParametersOP which is constructed from the given input
	/// stream (istream).

	utility::pointer::ReferenceCountCOP
	create_resource(
		basic::resource_manager::ResourceManager & resource_manager,
		utility::tag::TagCOP resource_tag,
		std::string const & input_id,
		std::istream & istream
	) const override;

	static
	SurfaceParametersOP
	create_surface_params(
		std::string const & input_id,
		std::istream & istream
	);

	static
	std::string
	classname();

	static
	void
	provide_xml_schema(
		utility::tag::XMLSchemaDefinition & xsd
	);

};

} // namespace surface_docking
} // namespace protocols

#endif //INCLUDED_protocols_surface_docking_SurfaceVectorLoader_HH
