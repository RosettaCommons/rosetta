// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceOptionsFactory.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceOptionsFactory_hh
#define INCLUDED_basic_resource_manager_ResourceOptionsFactory_hh

// package headers
#include <basic/resource_manager/ResourceOptionsCreator.fwd.hh>
#include <basic/resource_manager/ResourceOptions.fwd.hh>

//project headers

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//C++ headers
#include <map>

namespace basic {
namespace resource_manager {

/// @brief The %ResourceOptionsFactory class is responsible for maintaining the map
/// between the names of the ResourceOptions classes (strings) and the ResourceOptionsCreator
/// classes that are responsible for instantiating the ResourceOption classes.  This
/// is a singleton class.  It is initialized at load time -- or at least after the call to
/// devel::init( argc, argv ) -- and populated with the help of ResourceOptionsRegistrator
/// instances.
class ResourceOptionsFactory
{
public:
	virtual ~ResourceOptionsFactory();

	/// @brief Create an instance of a ResourceOptions class given the name of the class
	/// (options_type) and initialize it from the given tag object.  Throws an
	/// EXCN_Msg_Exception exception if the options_type string is not recognized.
	ResourceOptionsOP
	create_resource_options(
		std::string const & options_type,
		utility::tag::TagPtr tag
	) const;

	/// @brief Singleton accessor function.  Return the globally-unique instance of the class.
	static
	ResourceOptionsFactory *
	get_instance();

	/// @brief Register the given ResourceOptionsCreator object with the factory.  Invoked by an
	/// instance of the templated ResourceOptionsRegistrator class.
	void
	factory_register( ResourceOptionsCreatorOP creator );

	/// @brief Instruct the %ResourceOptionsFactory to throw an exception if two ResourceOptionCreators
	/// are registered with the factory that report the same "options_type" string instead of invoking
	/// utility_exit_with_message.
	void
	set_throw_on_double_registration();

private:

	/// @brief Singleton private constructor
	ResourceOptionsFactory();

private:
	static ResourceOptionsFactory * instance_;

	bool throw_on_double_registration_;

	typedef std::map< std::string, ResourceOptionsCreatorOP > ResourceOptionsCreatorMap;
	ResourceOptionsCreatorMap creator_map_;

};

} // namespace resource_manager
} // namespace basic



#endif //INCLUDED_basic_resource_manager_ResourceOptionsFactory_hh
