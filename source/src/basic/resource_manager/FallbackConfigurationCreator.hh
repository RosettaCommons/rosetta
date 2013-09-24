// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/FallbackConfigurationCreator.hh
/// @brief  Declaration of the FallbackConfigurationCreator class
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

/// @brief The %FallbackConfigurationCreator plays the role in the
/// ResourceManager framework of gluing a "resource_description" string
/// and a FallbackConfiguration together.

/// @details The ResourceManager looks to a FallbackConfiguration when a resource is
/// requested by a resource description string but no resource definition file
/// has been provided.  In such a case, the FallbackConfiguration will provide
/// the information the ResourceManager needs to create the resource.  It's the
/// %FallbackConfigurationCreator's job to inform the ResourceManager which
/// FallbackConfiguration to talk to.
///
/// Each class derived from the %FallbackConfigurationCreator will
/// instantiate a single FallbackConfiguration and act to pair a string,
/// a "resource description," with that FallbackConfiguration.  For example
/// "LoopFile" as a resource description will be paired by the
/// LoopFileFallbackConfigurationCreator with the LoopFileFallbackConfiguration.
/// Multiple resource descriptions can be paired with a single FallbackConfiguration.
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
