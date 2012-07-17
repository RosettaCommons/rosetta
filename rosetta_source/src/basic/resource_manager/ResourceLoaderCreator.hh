// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

namespace basic {
namespace resource_manager {

class ResourceLoaderCreator : public utility::pointer::ReferenceCount
{
public:
	virtual
	~ResourceLoaderCreator();

	virtual
	ResourceLoaderOP
	create_resource_loader() const = 0;

	virtual
	std::string loader_type() const = 0;

};

} // namespace resource_manager
} // namespace basic



#endif //INCLUDED_basic_resource_manager_ResourceLoaderCreator_hh
