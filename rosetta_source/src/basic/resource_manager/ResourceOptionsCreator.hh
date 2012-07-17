// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceOptionsCreator.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceOptionsCreator_HH
#define INCLUDED_basic_resource_manager_ResourceOptionsCreator_HH

// unit headers
#include <basic/resource_manager/ResourceOptionsCreator.fwd.hh>

// package headers
#include <basic/resource_manager/ResourceOptions.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace basic {
namespace resource_manager {

class ResourceOptionsCreator : public utility::pointer::ReferenceCount
{
public:
	virtual ~ResourceOptionsCreator();

	virtual std::string options_type() const = 0;
	virtual ResourceOptionsOP create_options() const = 0;

};

} // namespace resource_manager
} // namespace basic



#endif //INCLUDED_basic_resource_manager_ResourceOptionsCreator_hh
