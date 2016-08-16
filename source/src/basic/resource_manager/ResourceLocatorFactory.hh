// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceLocatorFactory.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceLocatorFactory_hh
#define INCLUDED_basic_resource_manager_ResourceLocatorFactory_hh


//package headers
#include <basic/resource_manager/ResourceLocatorCreator.fwd.hh>
#include <basic/resource_manager/ResourceLocator.fwd.hh>

//utility headers
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//C++ headers
#include <map>

namespace basic {
namespace resource_manager {

/// @brief The %ResourceLocatorFactory instantiates ResourceLocator objects given their corresponding
/// locator-type strings.  It uses the load-time factory registration scheme, meaning that it is a
/// singleton and takes an instance of a Creator object (a ResourceLocatorCreator) in its
/// "factory_register" method.  Templated instances of the ResourceLocatorRegistrator classes should
/// be placed in the library init.cc files  (e.g. core/init/init.cc or
/// protocols/init/init.ResourceLocatorRegistrators.ihh)
class ResourceLocatorFactory : public utility::SingletonBase< ResourceLocatorFactory >
{
public:
	friend class utility::SingletonBase< ResourceLocatorFactory >;

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~ResourceLocatorFactory();

	/// @brief Create a ResourceLocator given its locator_type, giving the newly created instance the name locator_tag
	ResourceLocatorOP
	create_resource_locator(
		std::string const & locator_type,
		std::string const & locator_tag,
		utility::tag::TagCOP tags
	) const;

	/// @brief This function is called on the singleton instance to give a ResourceLocatorCreator to the
	/// factory, usually through the constructor of a ResourceLocatorRegistrator class.
	void
	factory_register( ResourceLocatorCreatorOP creator );

	/// @brief Only useful for unit testing.  Since factory registration happens (sometimes) at
	/// load time, there may be no one to catch a thrown exception in the event of a name collision
	/// between two ResourceLocatorCreators that register for the same name.
	void
	set_throw_on_double_registration();

private:

	/// @brief Singleton has a private constructor
	ResourceLocatorFactory();

private:

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static ResourceLocatorFactory * create_singleton_instance();

	bool throw_on_double_registration_;
	std::map< std::string, ResourceLocatorCreatorOP > creator_map_;

};

} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_manager_ResourceLocatorFactory_hh
