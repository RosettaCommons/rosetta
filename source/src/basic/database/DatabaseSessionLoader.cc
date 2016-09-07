// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/database/DatabaseSessionLoader.cc
/// @brief  load the database session object.
/// @author Tim Jacobs

//unit headers
#include <basic/database/DatabaseSessionLoader.hh>
#include <basic/database/DatabaseSessionLoaderCreator.hh>
#include <basic/database/DatabaseSessionOptionsCreator.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/excn/Exceptions.hh>

//package headers
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>

//C++ headers
#include <istream>

namespace basic {
namespace database {

basic::resource_manager::ResourceOP
DatabaseSessionLoader::create_resource(
	basic::resource_manager::ResourceOptions const & options,
	basic::resource_manager::LocatorID const &,
	std::istream & /*istream*/
) const {
	if ( ! dynamic_cast< DatabaseSessionOptions const * > ( &options ) ) {
		throw utility::excn::EXCN_Msg_Exception( "DatabaseSessionLoader expected to be given a DatabaseSessionOptions object, " \
			"but was given a non-DatabaseSessionOptions object of type '" + options.type() + "', which has the name '" + options.name() + "'." );
	}
	DatabaseSessionOptions const & database_options = static_cast< DatabaseSessionOptions const & > ( options );
	return database_options.database_session();
}


///// DatabaseSessionOptionsCreator /////
DatabaseSessionOptionsCreator::DatabaseSessionOptionsCreator() {}

DatabaseSessionOptionsCreator::~DatabaseSessionOptionsCreator() = default;

basic::resource_manager::ResourceOptionsOP
DatabaseSessionOptionsCreator::create_options() const {
	return basic::resource_manager::ResourceOptionsOP( new DatabaseSessionOptions );
}

std::string
DatabaseSessionOptionsCreator::options_type() const {
	return "DatabaseSessionOptions";
}

//// DatabaseSessionLoaderCreator
basic::resource_manager::ResourceLoaderOP
DatabaseSessionLoaderCreator::create_resource_loader() const
{
	return basic::resource_manager::ResourceLoaderOP( new DatabaseSessionLoader() );
}

std::string DatabaseSessionLoaderCreator::loader_type() const
{
	return "DatabaseSession";
}


} // database
} // basic
