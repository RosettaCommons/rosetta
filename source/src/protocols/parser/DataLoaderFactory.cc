// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/parser/DataLoaderFactory.cc
/// @brief  Implementation of the factory class for the parser's DataLoader classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/parser/DataLoaderFactory.hh>

// Package Headers
#include <protocols/parser/DataLoader.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace protocols {
namespace parser {

DataLoaderFactory::~DataLoaderFactory() {}

void
DataLoaderFactory::factory_register( DataLoaderCreatorOP creator )
{
	//std::cout << "DataLoaderFactory::factory_register of " << creator->keyname() << std::endl;

	runtime_assert( creator != 0 );
	std::string const loader_type( creator->keyname() );
	if ( loader_type == "UNDEFINED NAME" ) {
		utility_exit_with_message("Can't map derived DataLoader with undefined type name.");
	}
	if ( dataloader_creator_map_.find( loader_type ) != dataloader_creator_map_.end() ) {
		utility_exit_with_message("DataLoaderFactory::factory_register already has a DataLoaderCreator with name \"" + loader_type + "\".  Conflicting Filter names" );
	}
	dataloader_creator_map_[ loader_type ] = creator;

}

/// @brief Create a DataLoader given its identifying string
DataLoaderOP
DataLoaderFactory::newDataLoader( std::string const & loader_type ) const
{

	//std::cout << "DataLoaderFactory::newDataLoader of " << loader_type << std::endl;

	LoaderMap::const_iterator iter( dataloader_creator_map_.find( loader_type ) );
	if ( iter != dataloader_creator_map_.end() ) {
		if ( ! iter->second ) {
			utility_exit_with_message( "Error: DataLoaderCreatorOP prototype for " + loader_type + " is NULL!" );
		}
		return iter->second->create_loader();
	} else {
		utility_exit_with_message( loader_type + " is not known to the DataLoaderFactory. Was it registered via a DataLoaderRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
		return NULL;
	}

}

DataLoaderFactory::LoaderMap const &
DataLoaderFactory::loader_map() const
{
	return dataloader_creator_map_;
}

std::string
DataLoaderFactory::data_loader_ct_namer( std::string const & loader_name )
{
	return "data_loader_" + loader_name + "_type";
}


DataLoaderFactory::DataLoaderFactory() {}

} //namespace parser
} //namespace protocols
