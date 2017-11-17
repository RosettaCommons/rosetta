// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/LoopHashLibraryLoader.cc
/// @brief Load the Loop Hash library using the resource manager
/// @author Tim Jacobs

//unit headers
#include <protocols/loophash/LoopHashLibraryLoader.hh>
#include <protocols/loophash/LoopHashLibraryLoaderCreator.hh>
#include <protocols/loophash/LoopHashLibraryOptionsCreator.hh>
#include <protocols/loophash/LoopHashLibrary.hh>

//package headers
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>

//C++ headers
#include <istream>

namespace protocols {
namespace loophash {

/// @details Return an owning pointer to a newly constructed default instance of LoopHashLibraryLoader.
basic::resource_manager::ResourceLoaderOP
LoopHashLibraryLoaderCreator::create_resource_loader() const
{
	return basic::resource_manager::ResourceLoaderOP( new LoopHashLibraryLoader() );
}

/// @details Return a string specifying the type of %ResourceLoader to create (LoopHashLibrary).
std::string LoopHashLibraryLoaderCreator::loader_type() const
{
	return "LoopHashLibrary";
}


LoopHashLibraryLoader::LoopHashLibraryLoader() {}

/// @details Ensure the %ResourceOptions is a LoopHashLibraryOptions instance and construct a new LoopHashLibrary from
/// it.  The locator_id and istream are not used.
/// @throws utility::excn::EXCN_Msg_Exception
basic::resource_manager::ResourceOP
LoopHashLibraryLoader::create_resource(
	basic::resource_manager::ResourceOptions const & options,
	basic::resource_manager::LocatorID const &,
	std::istream &
) const {
	if ( ! dynamic_cast< LoopHashLibraryOptions const * > ( &options ) ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "LoopHashLibraryLoader expected to be given a LoopHashLibraryOptions object, " \
			"but was given a non-LoopHashLibraryOptions object of type '" + options.type() + "', which has the name '" + options.name() + "'." );
	}
	LoopHashLibraryOptions const & lh_opts = static_cast< LoopHashLibraryOptions const & > ( options );
	LoopHashLibraryOP lh_library( new LoopHashLibrary( lh_opts.loop_sizes() ) );
	lh_library->load_mergeddb();
	return lh_library;
}

} // namespace
} // namespace
