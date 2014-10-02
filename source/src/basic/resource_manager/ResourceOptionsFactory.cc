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
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>

//C++ headers
#include <sstream>

// Singleton instance and mutex static data members
namespace utility {

using basic::resource_manager::ResourceOptionsFactory;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< ResourceOptionsFactory > ::singleton_mutex_;
template <> std::atomic< ResourceOptionsFactory * > utility::SingletonBase< ResourceOptionsFactory >::instance_( 0 );
#else
template <> ResourceOptionsFactory * utility::SingletonBase< ResourceOptionsFactory >::instance_( 0 );
#endif

}

namespace basic {
namespace resource_manager {

ResourceOptionsFactory *
ResourceOptionsFactory::create_singleton_instance()
{
	return new ResourceOptionsFactory;
}

ResourceOptionsFactory::~ResourceOptionsFactory() {}

ResourceOptionsOP
ResourceOptionsFactory::create_resource_options(
	std::string const & options_type,
	utility::tag::TagCOP tag
) const
{
	ResourceOptionsCreatorMap::const_iterator iter = creator_map_.find( options_type );
	if ( iter == creator_map_.end() ) {
		std::stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized ResourceOptions "
			<< "'" << options_type << "'." << std::endl
			<< "Check the spelling or "
			<< "register a new ResourceOptionsCreator with the ResourceOptionsFactory." << std::endl
			<< "Known ResourceOptions types are:" << std::endl;

		BOOST_FOREACH(const ResourceOptionsCreatorMap::value_type& type, creator_map_){
			error_msg << "\t" << type.first << std::endl;
		}

		throw utility::excn::EXCN_Msg_Exception(error_msg.str());
	}
	ResourceOptionsOP resource_options = (*iter).second->create_options();
	resource_options->parse_my_tag( tag );
	return resource_options;
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

/// @details Only useful for unit testing.  Since factory registration can happen any time between
/// load time and the call to devel::init, there may be no one to catch a thrown exception in the
/// event of a name collision between two ResourceOptionsCreators that register for the same name,
/// and so the default behavior is to call utility_exit_with_message.  For the sake of unit testing, however,
/// it is useful to change the default behavior so that exceptions can be thrown and then caught.
void
ResourceOptionsFactory::set_throw_on_double_registration() { throw_on_double_registration_ = true; }

ResourceOptionsFactory::ResourceOptionsFactory() :
	throw_on_double_registration_( false )
{}


} // namespace resource_manager
} // namespace basic
