// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/looops/SurfaceVectorFallbackConfiguration.hh
/// @author Michael Pacella mpacella88@gmail.com

#ifndef INCLUDED_protocols_surface_docking_surface_vector_fallback_configuration_HH
#define INCLUDED_protocols_surface_docking_surface_vector_fallback_configuration_HH

// Unit Header
#include <protocols/surface_docking/SurfaceVectorFallbackConfiguration.fwd.hh>
#include <basic/resource_manager/FallbackConfiguration.hh>
#include <basic/resource_manager/types.hh>

namespace protocols {
namespace surface_docking {

class SurfaceVectorFallbackConfiguration : public basic::resource_manager::FallbackConfiguration {
public:
	typedef basic::resource_manager::ResourceDescription ResourceDescription;
public:
	SurfaceVectorFallbackConfiguration();

	virtual
	bool
	fallback_specified( ResourceDescription const & desc ) const;

	virtual
	basic::resource_manager::LoaderType
	get_resource_loader( ResourceDescription const & desc ) const;

	virtual
	basic::resource_manager::LocatorID
	get_locator_id( ResourceDescription const & desc ) const;

	virtual
	basic::resource_manager::ResourceOptionsOP
	get_resource_options( ResourceDescription const & desc ) const;

	virtual
	std::string
	could_not_create_resource_error_message( ResourceDescription const & desc ) const;

private:

	basic::resource_manager::LocatorID get_surface_vector_filename_from_options() const;

};

} // namespace surface_docking
} // namespace protocols

#endif // INCLUDED_protocols_surface_docking_surface_vector_fallback_configuration_HH
