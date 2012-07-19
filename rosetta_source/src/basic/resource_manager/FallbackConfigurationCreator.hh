// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/FallbackConfigurationCreator.hh
/// @brief  
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_basic_resource_manager_fallback_configuration_creator_HH
#define INCLUDED_basic_resource_manager_fallback_configuration_creator_HH

//unit headers
#include <basic/resource_manager/FallbackConfigurationCreator.fwd.hh>

// package headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <string>

namespace basic {
namespace resource_manager {

class FallbackConfigurationCreator : public utility::pointer::ReferenceCount
{
public:
	virtual
	~FallbackConfigurationCreator();

	virtual
	FallbackConfigurationOP
	create_fallback_configuration() const = 0;

	virtual
	std::string resource_description() const = 0;

};


} // namespace resource_manager
} // namespace basic

#endif // INCLUDED_basic_resource_manager_fallback_configuration_creator_HH
