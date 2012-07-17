// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceLoader.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceLoader_hh
#define INCLUDED_basic_resource_manager_ResourceLoader_hh

//unit headers
#include <basic/resource_manager/ResourceLoader.fwd.hh>

//package headers
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <istream>

namespace basic {
namespace resource_manager {

class ResourceLoader : public utility::pointer::ReferenceCount
{
public:
	virtual ~ResourceLoader();

	/// @brief Create a resource, held in an owning pointer, of any type
	/// which will be stored and whose lifetime will be governed by the\
	/// ResourceManager
	virtual
	ResourceOP
	create_resource(
		ResourceOptions const & options,
		LocatorID const & locator_id,
		std::istream & istream
	) const = 0;

	virtual
	ResourceOptionsOP
	default_options() const = 0;
};

} // namespace resource_manager
} // namespace basic



#endif //INCLUDED_basic_resource_manager_ResourceLoader_hh
