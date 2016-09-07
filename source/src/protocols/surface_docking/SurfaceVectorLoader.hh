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

//project headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

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
	
	utility::pointer::ReferenceCountOP
	create_resource(
		basic::resource_manager::ResourceOptions const & options,
		basic::resource_manager::LocatorID const & locator_id,
		std::istream & istream
	) const override;

	/// @brief Returns the default options for SurfaceParameters
	
	basic::resource_manager::ResourceOptionsOP
	default_options() const override;

};

} // namespace surface_docking
} // namespace protocols

#endif //INCLUDED_protocols_surface_docking_SurfaceVectorLoader_HH
