// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

utility::pointer::ReferenceCountOP
LoopsFileLoader::create_resource(
	basic::resource_manager::ResourceOptions const & options,
	basic::resource_manager::LocatorID const & locator_id,
	std::istream & istream
) const
{
	LoopsFileOptions const * loops_opts_ptr = dynamic_cast< LoopsFileOptions const * > ( &options );
	if ( ! loops_opts_ptr ) {
		throw utility::excn::EXCN_Msg_Exception( "LoopsFileLoader expected to be given a LoopsFileOptions object, " \
			"but was given a non-LoopsFileOptions object of type '" + options.type() + "', which has the name '" + options.name() + "'." );
	}
	LoopsFileIO lfio;
	LoopsFileDataOP lfd = lfio.read_loop_file_stream( istream, locator_id, loops_opts_ptr->prohibit_single_residue_loops() );
	return lfd;
}

basic::resource_manager::ResourceOptionsOP
LoopsFileLoader::default_options() const
{
	return new LoopsFileOptions;
}

basic::resource_manager::ResourceLoaderOP LoopsFileLoaderCreator::create_resource_loader() const
{
	return new LoopsFileLoader();
}

std::string LoopsFileLoaderCreator::loader_type() const
{
	return "LoopsFile";
}

} // namespace loops
} // namespace protocols
