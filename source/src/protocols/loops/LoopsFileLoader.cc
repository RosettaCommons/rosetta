// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/LoopsFileLoader.cc
/// @brief
/// @author

//unit headers
#include <protocols/loops/LoopsFileLoader.hh>
#include <protocols/loops/LoopsFileLoaderCreator.hh>

//package headers
#include <protocols/loops/LoopsFileOptions.hh>
#include <protocols/loops/LoopsFileIO.hh>

//utility headers
#include <utility/excn/Exceptions.hh>

// numeric headers

//C++ headers
#include <istream>

namespace protocols {
namespace loops {

LoopsFileLoader::LoopsFileLoader() {}
LoopsFileLoader::~LoopsFileLoader() {}

/// @details Ensure the %ResourceOptions is a LoopsFileOptions instance and construct a new LoopsFileData from the
/// istream and the options.  The locator_id is used solely for reporting accurate error messages.
/// @throws utility::excn::EXCN_Msg_Exception
utility::pointer::ReferenceCountOP
LoopsFileLoader::create_resource(
	basic::resource_manager::ResourceOptions const & options,
	basic::resource_manager::LocatorID const & locator_id,
	std::istream & istream
) const
{
	if ( ! dynamic_cast< LoopsFileOptions const * > ( &options ) ) {
		throw utility::excn::EXCN_Msg_Exception( "LoopsFileLoader expected to be given a LoopsFileOptions object, " \
			"but was given a non-LoopsFileOptions object of type '" + options.type() + "', which has the name '" + options.name() + "'." );
	}
	LoopsFileOptions const & loops_opts = static_cast< LoopsFileOptions const & > ( options );
	LoopsFileIO lfio;
	LoopsFileDataOP lfd = lfio.read_loop_file_stream( istream, locator_id, loops_opts.prohibit_single_residue_loops() );
	return lfd;
}

/// @details Return an owning pointer to a newly constructed default instance of LoopsFileOptions.
basic::resource_manager::ResourceOptionsOP
LoopsFileLoader::default_options() const
{
	return basic::resource_manager::ResourceOptionsOP( new LoopsFileOptions );
}

/// @details Return an owning pointer to a newly constructed default instance of LoopsFileLoader.
basic::resource_manager::ResourceLoaderOP LoopsFileLoaderCreator::create_resource_loader() const
{
	return basic::resource_manager::ResourceLoaderOP( new LoopsFileLoader() );
}

/// @details Return a string specifying the type of %ResourceLoader to create (LoopsFile).
std::string LoopsFileLoaderCreator::loader_type() const
{
	return "LoopsFile";
}

} // namespace loops
} // namespace protocols
