// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceManager.cc
/// @brief
/// @author


//unit headers
#include <basic/resource_manager/ResourceManager.hh>

//project headres
#include <basic/resource_manager/ResourceManagerFactory.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/thread/threadsafe_creation.hh>

//C++ headers
#include <sstream>
#include <iomanip>


// // Singleton instance and mutex static data members
// namespace utility {
//
// using basic::resource_manager::ResourceManager;
//
// #if defined MULTI_THREADED && defined CXX11
// template <> std::mutex utility::SingletonBase< ResourceManager > ::singleton_mutex_;
// template <> std::atomic< ResourceManager * > utility::SingletonBase< ResourceManager >::instance_( 0 );
// #else
// template <> ResourceManager * utility::SingletonBase< ResourceManager >::instance_( 0 );
// #endif
//
// }

namespace basic {
namespace resource_manager {

using std::stringstream;
using std::endl;
using std::setw;

ResourceManager *
ResourceManager::get_instance()
{
	static ResourceManager * instance_( 0 );
	// NOT THREAD SAFE!
	if ( ! instance_ ) {
		instance_ = ResourceManagerFactory::get_instance()->create_resource_manager_from_options_system();
	}
	return instance_;
}

//ResourceManager *
//ResourceManager::create_singleton_instance()
//{
// return ResourceManagerFactory::get_instance()->create_resource_manager_from_options_system();
//}

ResourceManager::ResourceManager() {}

ResourceManager::~ResourceManager() {}


void
ResourceManager::add_resource(
	ResourceTag const & resource_tag,
	ResourceOP resource
) {
	resources_[resource_tag] = resource;
}

bool
ResourceManager::has_resource(
	ResourceTag const & resource_tag
) const {
	ResourcesMap::const_iterator resource(
		resources_.find(resource_tag));
	return resource != resources_.end();
}

ResourceOP
ResourceManager::find_resource(
	ResourceTag const & resource_tag
) {
	ResourcesMap::const_iterator resource(
		resources_.find(resource_tag));
	if ( resource == resources_.end() ) {
		stringstream err_msg;
		err_msg
			<< "Unable to find resource for the resource tag '"
			<< resource_tag << "'" << endl;
		utility_exit_with_message(err_msg.str());
	}
	return resource->second;
}

void
ResourceManager::clear()
{
	resources_.clear();
}


void
ResourceManager::free_resource(
	ResourceTag const & resource_tag
) {
	if ( !has_resource(resource_tag) ) {
		stringstream err_msg;
		err_msg
			<< "Attempting to free resource with tag '" << resource_tag << "' "
			<< "but the resource does not exist.";
		utility_exit_with_message(err_msg.str());
	} else {
		resources_.erase(resource_tag);
	}
}

void
ResourceManager::show(
	std::ostream & out
) const {
	out
		<< "ResourceManager.resources:" << endl
		<< std::setiosflags(std::ios::left) << setw(16) << "ResourceTag" << "ResourceExists" << endl;
	for (
			ResourceManager::ResourcesMap::const_iterator
			r = resources_.begin(), re = resources_.end(); r != re; ++r ) {
		out
			<< std::setiosflags(std::ios::left) << setw(16) << r->first
			<< (r->second ? "true" : "false") << endl;
	}
}

std::ostream &
operator<<(
	std::ostream & out,
	ResourceManager const & resource_manager
) {
	resource_manager.show(out);
	return out;
}


} // namespace resource_manager
} // namespace basic
