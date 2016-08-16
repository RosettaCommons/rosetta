// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

/// @brief The ResourceLoader is responsible for instantiating a Resource object
/// and initializing it.  In order to do so, the ResourceLoader is given an input
/// stream and a ResourceOptions object.  Note that the ResourceOptions object has
/// to be of the right type, or the ResourceLoader will not be able to read the
/// data that it needs out of it.  If the ResourceLoader is given the wrong kind
/// of ResourceOptions object, it will throw an exception.
class ResourceLoader : public utility::pointer::ReferenceCount
{
public:
	virtual ~ResourceLoader();

	/// @brief Create a resource, held in an owning pointer, of any type
	/// which will be stored and whose lifetime will be governed by the
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
