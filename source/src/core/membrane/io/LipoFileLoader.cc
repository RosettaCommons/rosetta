// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileLoader.cc
///
/// @brief      Loader class - load object that stores Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_LipoFileLoader_cc
#define INCLUDED_core_membrane_io_LipoFileLoader_cc

// Unit Headers
#include <core/membrane/io/LipoFileLoader.hh>
#include <core/membrane/io/LipoFileLoaderCreator.hh>

// Project Headers
#include <core/conformation/membrane/LipidAccInfo.hh>
#include <core/conformation/membrane/Exceptions.hh>

#include <core/membrane/io/LipoFileOptions.hh>

// Utility Headers
#include <utility/excn/Exceptions.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace membrane {
namespace io {

/// @brief Constructor
LipoFileLoader::LipoFileLoader() {}
    
/// @brief Destructor
LipoFileLoader::~LipoFileLoader() {}

/// @brief Create a Lipid Builder Object from Data File
utility::pointer::ReferenceCountOP
LipoFileLoader::create_resource(
	basic::resource_manager::ResourceOptions const & options,
	basic::resource_manager::LocatorID const & locator_id,
	std::istream &
	) const
{
	using namespace core::conformation::membrane;
	using namespace core::membrane::io;

	// Create and initialize a lipid object
	LipidAccInfoOP lips_exp = new LipidAccInfo( locator_id );

	// Cast generic options to lipofile opts
	if ( ! dynamic_cast< LipoFileOptions const * > ( &options ) ) {
		throw utility::excn::EXCN_Msg_Exception("LipoFileLoader excpected to be given a LipoFileOptions but was given "
				"a non-LipoFileOptions object of type " + options.type() + "' which has the name '"
				+ options.name() + "'.");
	}

	return lips_exp;
}
    
/// @brief Import a Default Options Object
basic::resource_manager::ResourceOptionsOP
LipoFileLoader::default_options() const
{
	return new LipoFileOptions;
}

/// @brief Return new resource loader for lipo file
basic::resource_manager::ResourceLoaderOP
LipoFileLoaderCreator::create_resource_loader() const
{
	return new LipoFileLoader();
}

/// @brief Lipo File Loader Type
std::string LipoFileLoaderCreator::loader_type() const
{
	return "LipoFile";
}

} // io
} // memrbane
} // core

#endif // INCLUDED_core_membrane_io_LipoFileLoader_cc


