// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceLoaderFactory.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceLoaderFactory_hh
#define INCLUDED_basic_resource_manager_ResourceLoaderFactory_hh

// Unit headers
#include <basic/resource_manager/ResourceLoaderCreator.fwd.hh>
#include <basic/resource_manager/ResourceLoader.fwd.hh>

// Package headers
#include <basic/resource_manager/ResourceManager.fwd.hh>

//utility headers
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

//C++ headers
#include <list>
#include <map>
#include <string>

#ifdef MULTI_THREADED
#include <atomic>
#include <mutex>
#endif

namespace basic {
namespace resource_manager {

/// @brief Instantiates ResourceLoaders.  Creators may be registered with
/// the Factory at any point, though it is recommended they be registered
/// at load time.  If two Creators are registered and they both give the
/// same name for the ResourceLoader they say they will instantiate, then
/// the Factory will exit with an error message.
class ResourceLoaderFactory : public utility::SingletonBase< ResourceLoaderFactory >
{
public:
	friend class utility::SingletonBase< ResourceLoaderFactory >;
public:
	/// @brief Instantiates a resource loader ofa  given type; throws
	/// an exception if no loader with this type has been previously
	/// registered.
	ResourceLoaderOP
	create_resource_loader(
		std::string const & loader_type
	) const;

	/// @brief Returns true if a resource loader of the given type has been
	/// registered with the factory
	bool
	has_resource_loader(
		std::string const & loader_type
	) const;

	/// @brief Return a list of all the resource loaders available
	std::list< std::string >
	available_resource_loaders() const;

	/// @brief Register a ResourceLoaderCreator with the factory.  The factory
	/// asks the Creator for the name of the ResourceLoader that it will create;
	/// if another Creator has already registered with the factory that proports
	/// to instantiate another ResourceLoader with the same name, then it will
	/// exit with an error message, or, if set_throw_on_double_registration()
	/// has previously been called, throw an exception.
	void
	factory_register( ResourceLoaderCreatorOP creator );

	/// @brief Only useful for unit testing.  Since factory registration happens (sometimes) at
	/// load time, there may be no one to catch a thrown exception in the event of a name collision
	void
	set_throw_on_double_registration();

	/// @brief Read access to the set of creators that the factry holds; used to create the
	/// schema definition for the ResourceManager.
	std::map< std::string, ResourceLoaderCreatorOP > const &
	loader_map() const;

	/// @brief The name mangling function used by ResourceLoaders when they define their complexTypes
	/// in the XML Schema.
	static
	std::string
	complex_type_name_for_loader( std::string const & loader_name );

private:

	/// singleton has a private constructor
	ResourceLoaderFactory();

	ResourceLoaderFactory( ResourceLoaderFactory const & ) = delete;
	ResourceLoaderFactory & operator=( ResourceLoaderFactory const & ) = delete;

private:

	bool throw_on_double_registration_;
	std::map< std::string, ResourceLoaderCreatorOP > creator_map_;

};

} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_manager_ResourceLoaderFactory_hh
