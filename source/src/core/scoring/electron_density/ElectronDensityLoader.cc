// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/electron_density/ElectronDensityLoader.cc
/// @brief  Options for constructing an electron density map
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <core/scoring/electron_density/ElectronDensityLoader.hh>
#include <core/scoring/electron_density/ElectronDensityLoaderCreator.hh>
#include <basic/resource_manager/ResourceOptions.hh>

// Project Headers
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/ElectronDensityOptions.hh>

// Platform Headers
#include <core/types.hh>
#include <utility/excn/Exceptions.hh>

// C++ Headers
#include <string>

namespace core {
namespace scoring {
namespace electron_density {

using basic::resource_manager::ResourceOP;
using basic::resource_manager::ResourceLoaderOP;
using basic::resource_manager::ResourceOptions;
using basic::resource_manager::ResourceOptionsOP;
using basic::resource_manager::LocatorID;
using std::istream;
using std::string;

///// ElectronDensityLoaderCreator /////
ElectronDensityLoaderCreator::ElectronDensityLoaderCreator() {}

ElectronDensityLoaderCreator::~ElectronDensityLoaderCreator() {}

ResourceLoaderOP
ElectronDensityLoaderCreator::create_resource_loader() const {
	return ResourceLoaderOP( new ElectronDensityLoader );
}

string
ElectronDensityLoaderCreator::loader_type() const {
	return "ElectronDensity";
}

//// ElectronDensityLoader /////
ElectronDensityLoader::ElectronDensityLoader() {}

ElectronDensityLoader::~ElectronDensityLoader() {}

ElectronDensityLoader::ElectronDensityLoader(
	ElectronDensityLoader const &) : basic::resource_manager::ResourceLoader() {}

ResourceOP
ElectronDensityLoader::create_resource(
	ResourceOptions const & options,
	LocatorID const & locator_id,
	istream & istream
) const {

	if ( ! dynamic_cast< ElectronDensityOptions const * >( &options ) ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,
			"ElectronDensityLoader expected to get a ElectronDensityOptions object, "
			"but was given a ResourceOptions of type '" + options.type() + "', "
			"which has the name '" + options.name() + "'." );
	}
	ElectronDensityOptions const & resource_options(
		static_cast< ElectronDensityOptions const & >( options ));

	ElectronDensityOP electron_density( new ElectronDensity() );

	electron_density->readMRCandResize(
		istream,
		locator_id,
		resource_options.get_mapreso(),
		resource_options.get_grid_spacing());

	return electron_density;
}

ResourceOptionsOP
ElectronDensityLoader::default_options(
) const {
	return ResourceOptionsOP( new ElectronDensityOptions() );
}

} // namespace
} // namespace
} // namespace

