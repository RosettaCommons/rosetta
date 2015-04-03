// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceLocatorCreator.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceLocatorCreator_hh
#define INCLUDED_basic_resource_manager_ResourceLocatorCreator_hh

//unit headers
#include <basic/resource_manager/ResourceLocatorCreator.fwd.hh>

//package headers
#include <basic/resource_manager/ResourceLocator.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#ifdef WIN32
#include <string>
#endif

namespace basic {
namespace resource_manager {

/// @brief The %ResourceLocatorCreator class serves to link the name of
/// a locator type and the (derived) ResourceLocator class that's responsible
/// for retrieving data from a data store.
class ResourceLocatorCreator : public utility::pointer::ReferenceCount
{
public:
	virtual
	~ResourceLocatorCreator();

	virtual
	ResourceLocatorOP
	create_resource_locator() const = 0;

	virtual
	std::string locator_type() const = 0;

};

} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_manager_ResourceLocatorCreator_hh
