// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       SpanFileLoader.cc
///
/// @brief      Loader for SpanningTopology object - generating per-chain membrane spanning data
/// @details    Generates membrane spanning topology data from Octopus topology
///             prediction information. Topology can be generated at http://octopus.cbr.su.se/
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_SpanFileLoader_cc
#define INCLUDED_core_membrane_io_SpanFileLoader_cc

// Unit Headers
#include <core/membrane/io/SpanFileLoader.hh>
#include <core/membrane/io/SpanFileLoaderCreator.hh>

// Project Headers
#include <core/membrane/io/SpanFileOptions.hh>

#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Utility Headers
#include <utility/excn/Exceptions.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static thread_local basic::Tracer TR( "Loader debug tracer: " );

namespace core {
namespace membrane {
namespace io {

/// @brief Constructor
SpanFileLoader::SpanFileLoader() {}
    
/// @brief Destructor
SpanFileLoader::~SpanFileLoader() {}

/// @brief Create a topology object from spanfile
utility::pointer::ReferenceCountOP
SpanFileLoader::create_resource(
	basic::resource_manager::ResourceOptions const & options,
	basic::resource_manager::LocatorID const & locator_id,
	std::istream &
	) const
{

	using namespace core::membrane::io;
	using namespace core::conformation::membrane;

	// Load a topology object from spanfile
	SpanningTopologyOP topology = new SpanningTopology( locator_id );

	// Adding opts casting
	if ( ! dynamic_cast< SpanFileOptions const * > ( &options ) ) {
		throw utility::excn::EXCN_Msg_Exception("SpanFileLoader excpected to be given a SpanFileOptions but was given "
				"a non-SpanFileOptions object of type " + options.type() + "' which has the name '"
				+ options.name() + "'.");
	}

	// Return a newly initialized topology object
	return topology;
}

/// @brief Import Default Options Object
basic::resource_manager::ResourceOptionsOP
SpanFileLoader::default_options() const
{
	return new SpanFileOptions;
}

/// @brief Return new resource loader for span file
basic::resource_manager::ResourceLoaderOP
SpanFileLoaderCreator::create_resource_loader() const
{
	return new SpanFileLoader();
}

/// @brief Span File Loader type
std::string SpanFileLoaderCreator::loader_type() const
{
	return "SpanFile";
}

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_SpanFileLoader_cc

