// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffractionLoader.cc
/// @brief  Options for constructing fiber diffraction layer lines
/// @author Wojciech Potrzebowski and Ingemar Andre

// Unit Headers
#include <core/scoring/fiber_diffraction/FiberDiffractionLoader.hh>
#include <core/scoring/fiber_diffraction/FiberDiffractionLoaderCreator.hh>
#include <basic/resource_manager/ResourceOptions.hh>

// Project Headers
#include <core/scoring/fiber_diffraction/FiberDiffraction.hh>
#include <core/scoring/fiber_diffraction/FiberDiffractionOptions.hh>

// Platform Headers
#include <core/types.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <string>

namespace core {
namespace scoring {
namespace fiber_diffraction {

using basic::resource_manager::ResourceOP;
using basic::resource_manager::ResourceLoaderOP;
using basic::resource_manager::ResourceOptions;
using basic::resource_manager::ResourceOptionsOP;
using basic::resource_manager::LocatorID;
using std::istream;
using std::string;

FiberDiffractionLoaderCreator::FiberDiffractionLoaderCreator() = default;

FiberDiffractionLoaderCreator::~FiberDiffractionLoaderCreator() = default;

ResourceLoaderOP
FiberDiffractionLoaderCreator::create_resource_loader() const {
	return ResourceLoaderOP( new FiberDiffractionLoader );
}

string
FiberDiffractionLoaderCreator::loader_type() const {
	return "FiberDiffraction";
}

FiberDiffractionLoader::FiberDiffractionLoader() = default;

FiberDiffractionLoader::~FiberDiffractionLoader() = default;

FiberDiffractionLoader::FiberDiffractionLoader(
	FiberDiffractionLoader const &) : basic::resource_manager::ResourceLoader() {}

ResourceOP
FiberDiffractionLoader::create_resource(
	ResourceOptions const & options,
	LocatorID const & locator_id,
	istream & istream
) const {

	if ( ! dynamic_cast< FiberDiffractionOptions const * >( &options ) ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,
			"FiberDiffractionLoader expected to get a FiberDiffractionOptions object, "
			"but was given a ResourceOptions of type '" + options.type() + "', "
			"which has the name '" + options.name() + "'." );
	}

	auto const & resource_options(
		static_cast< FiberDiffractionOptions const & >( options ));

	FiberDiffractionOP fiber_diffraction( new FiberDiffraction() );

	fiber_diffraction->loadFiberDiffractionData(
		istream,
		locator_id,
		resource_options.get_c_repeat(),
		resource_options.get_res_high(),
		resource_options.get_res_low());

	return fiber_diffraction;
}

ResourceOptionsOP
FiberDiffractionLoader::default_options() const {
	return ResourceOptionsOP( new FiberDiffractionOptions() );
}

} // namespace
} // namespace
} // namespace

