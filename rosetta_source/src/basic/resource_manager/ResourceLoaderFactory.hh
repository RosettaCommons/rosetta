// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceLoaderFactory.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceLoaderFactory_hh
#define INCLUDED_basic_resource_manager_ResourceLoaderFactory_hh


//package headers
#include <basic/resource_manager/ResourceLoaderCreator.fwd.hh>
#include <basic/resource_manager/ResourceLoader.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <list>
#include <map>
#include <string>

namespace basic {
namespace resource_manager {

class ResourceLoaderFactory : public utility::pointer::ReferenceCount
{
public:
	ResourceLoaderOP
	create_resource_loader(
		std::string const & loader_type
	) const;

	bool
	has_resource_loader(
		std::string const & loader_type
	) const;

	std::list< std::string >
	available_resource_loaders() const;

	static
	ResourceLoaderFactory *
	get_instance();

	void
	factory_register( ResourceLoaderCreatorOP creator );

private:

	/// singleton has a private constructor
	ResourceLoaderFactory();

private:
	static ResourceLoaderFactory * instance_;

	std::map< std::string, ResourceLoaderCreatorOP > creator_map_;

};

} // namespace resource_manager
} // namespace basic



#endif //INCLUDED_basic_resource_manager_ResourceLoaderFactory_hh
