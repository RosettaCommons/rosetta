// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceOptionsFactory.cc
/// @brief
/// @author

//unit headers
#include <basic/resource_manager/ResourceOptionsFactory.hh>

// package headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>
#include <basic/resource_manager/ResourceOptions.hh>

//project headers

//utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

//C++ headers

namespace basic {
namespace resource_manager {

ResourceOptionsFactory * ResourceOptionsFactory::instance_( 0 );

ResourceOptionsFactory::~ResourceOptionsFactory() {}

ResourceOptionsOP
ResourceOptionsFactory::create_resource_options(
	std::string const & options_type,
	utility::tag::TagPtr tag
) const
{
	std::map< std::string, ResourceOptionsCreatorOP >::const_iterator iter = creator_map_.find( options_type );
	if ( iter == creator_map_.end() ) {
		throw utility::excn::EXCN_Msg_Exception( "No ResourceOptionsCreator responsible for instantiating the ResourceOptions named " + options_type + " was found in the ResourceOptionsFactory.  Was it correctly registered?" );
	}
	ResourceOptionsOP resource_options = (*iter).second->create_options();
	resource_options->parse_my_tag( tag );
	return resource_options;
}

ResourceOptionsFactory *
ResourceOptionsFactory::get_instance()
{
	if ( ! instance_ ) {
		instance_ = new ResourceOptionsFactory;
	}
	return instance_;
}

void
ResourceOptionsFactory::factory_register( ResourceOptionsCreatorOP creator )
{
	std::string options_type = creator->options_type();
	std::map< std::string, ResourceOptionsCreatorOP >::const_iterator iter = creator_map_.find( options_type );
	if ( iter != creator_map_.end() ) {
		std::string errmsg( "Double registration of a ResourceOptionsCreator in the ResourceOptionsFactory, named " + options_type + ". Are there two registrators for this options object, or have you chosen a previously assigned name to a new resource option?" );
		if ( throw_on_double_registration_ ) {
			throw utility::excn::EXCN_Msg_Exception( errmsg );
		} else {
			utility_exit_with_message( errmsg );
		}
	}
	creator_map_[ options_type ] = creator;
}

void
ResourceOptionsFactory::set_throw_on_double_registration() { throw_on_double_registration_ = true; }

ResourceOptionsFactory::ResourceOptionsFactory() :
	throw_on_double_registration_( false )
{}


} // namespace resource_manager
} // namespace basic
