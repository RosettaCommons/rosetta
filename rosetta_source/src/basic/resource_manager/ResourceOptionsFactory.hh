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

class ResourceOptionsFactory : public utility::pointer::ReferenceCount
{
public:
	~ResourceOptionsFactory();

	ResourceOptionsOP
	create_resource_options(
		std::string const & options_type,
		utility::tag::TagPtr tag
	) const;

	static
	ResourceOptionsFactory *
	get_instance();

	void
	factory_register( ResourceOptionsCreatorOP creator );

	/// @brief Only useful for unit testing.  Since factory registration happens (sometimes) at
	/// load time, there may be no one to catch a thrown exception in the event of a name collision
	/// between two ResourceOptionsCreators that register for the same name. 
	void
	set_throw_on_double_registration();

private:

	/// singleton has a private constructor
	ResourceOptionsFactory();

private:
	static ResourceOptionsFactory * instance_;

	bool throw_on_double_registration_;
	std::map< std::string, ResourceOptionsCreatorOP > creator_map_;

};

} // namespace resource_manager
} // namespace basic



#endif //INCLUDED_basic_resource_manager_ResourceOptionsFactory_hh
