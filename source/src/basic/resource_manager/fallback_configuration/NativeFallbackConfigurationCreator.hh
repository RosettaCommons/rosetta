// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/fallback_configurations/NativeFallbackConfigurationCreator.hh
/// @author Matthew O'Meara mattjomeara@gmail.com

#ifndef INCLUDED_basic_resource_manager_fallback_configuration_NativeFallbackConfigurationCreator_HH
#define INCLUDED_basic_resource_manager_fallback_configuration_NativeFallbackConfigurationCreator_HH

//unit headers
#include <basic/resource_manager/FallbackConfigurationCreator.hh>

// package headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>

//C++ headers
#include <string>

namespace basic {
namespace resource_manager {
namespace fallback_configuration {

class NativeFallbackConfigurationCreator : public FallbackConfigurationCreator
{
public:
	virtual
	FallbackConfigurationOP
	create_fallback_configuration() const;

	virtual
	std::string resource_description() const;

};


} // namespace
} // namespace
} // namespace

#endif // include guard
