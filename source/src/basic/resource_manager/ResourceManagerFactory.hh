// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceManagerFactory.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceManagerFactory_HH
#define INCLUDED_basic_resource_manager_ResourceManagerFactory_HH

//unit headers
#include <basic/resource_manager/ResourceManagerFactory.fwd.hh>

// package headers
#include <basic/resource_manager/ResourceManagerCreator.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/factory/WidgetRegistrator.hh>

//C++ headers
#include <map>
#include <string>

namespace basic {
namespace resource_manager {

/// @brief A factory class for managing the instantiation of the singleton
/// ResourceManager: only one of the various derived classes will be instantiated.
/// Currently, it asks for the JD2ResourceManager; in the future, this should be
/// fixed so that it reads from the options system to figure out which ResourceManager
/// to instantiate.
class ResourceManagerFactory : public utility::SingletonBase< ResourceManagerFactory >
{
public:
	friend class utility::SingletonBase< ResourceManagerFactory >;

private:
	ResourceManagerFactory(); // singleton, private constructor

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static ResourceManagerFactory * create_singleton_instance();

public:

	/// Should only be called by the ResourceManager in its singleton construction!
	ResourceManager *
	create_resource_manager_from_options_system() const;

	void
	factory_register( ResourceManagerCreatorOP creator );

private:

	typedef std::map< std::string, ResourceManagerCreatorOP > ResourceManagerCreatorsMap;
	ResourceManagerCreatorsMap creators_map_;

};

/// @brief The %ResourceManagerRegistrator creates an instantiate of the templated ResourceManagerCreator
/// class and gives it to the ResourceManagerFactory.  A single instance for each ResourceManagerCreator
/// should be put in the appropriate init.cc file.
template < class T >
class ResourceManagerRegistrator : public utility::factory::WidgetRegistrator< ResourceManagerFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< ResourceManagerFactory, T > parent;
public:
	ResourceManagerRegistrator() : parent() {}
};


} // namespace resource_manager
} // namespace basic

#endif //INCLUDED_basic_resource_manager_ResourceManager_HH
