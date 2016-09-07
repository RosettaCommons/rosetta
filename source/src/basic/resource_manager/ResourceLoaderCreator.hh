// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceLoaderCreator.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceLoaderCreator_hh
#define INCLUDED_basic_resource_manager_ResourceLoaderCreator_hh

//unit headers
#include <basic/resource_manager/ResourceLoaderCreator.fwd.hh>

//package headers
#include <basic/resource_manager/ResourceLoader.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#ifdef WIN32
#include <string>
#endif

namespace basic {
namespace resource_manager {

/// @brief Instantiates a ResourceLoader as part of the ResourceLoaderFactory scheme.
/// Derived classes should be registered with the ResourceLoaderFactory in one of
/// the library init.cc files with a ResourceLoaderRegistrator
class ResourceLoaderCreator : public utility::pointer::ReferenceCount
{
public:

	~ResourceLoaderCreator() override;

	/// @brief Instantiate a ResourceLoader
	virtual
	ResourceLoaderOP
	create_resource_loader() const = 0;

	/// @brief Give the name of the ResourceLoader that this Creator will instantiate.
	virtual
	std::string loader_type() const = 0;

};

} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_manager_ResourceLoaderCreator_hh
