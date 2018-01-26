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
#include <basic/resource_manager/ResourceManager.fwd.hh>
#include <basic/resource_manager/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

//C++ headers
#include <istream>

namespace basic {
namespace resource_manager {

/// @brief The ResourceLoader is responsible for instantiating a Resource object
/// and initializing it.  In order to do so, the ResourceLoader is given an input
/// stream and a XML "Tag" object. The ResourceLoader can also request other
/// Resources from the ResourceManager in trying to construct a particular resource.
class ResourceLoader : public utility::pointer::ReferenceCount
{
public:
	virtual ~ResourceLoader();

	/// @brief Create a resource, held in an owning pointer, of any type
	/// which will be stored in the ResourceManager; a resource may depend
	/// on another resource held / managed by the ResourceManager. If a
	/// second (or third, or fourth, etc.) resource is requested during the
	/// construction of this resource, the resource manager will take note of
	/// their interdependency and preserve the second resource for as long
	/// as the first one has not been deallocated.
	virtual
	ResourceCOP
	create_resource(
		ResourceManager & resource_manager,
		utility::tag::TagCOP resource_tag,
		std::string const & input_id,
		std::istream & istream
	) const = 0;

};

} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_manager_ResourceLoader_hh
