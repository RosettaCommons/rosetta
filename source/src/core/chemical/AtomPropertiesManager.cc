// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/AtomPropertiesManager.cc
/// @brief   Method definitions for AtomPropertiesManager.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit headers
#include <core/chemical/AtomPropertiesManager.hh>
#include <core/chemical/AtomProperty.hh>


namespace core {
namespace chemical {

// Public methods /////////////////////////////////////////////////////////////
AtomProperty const &
AtomPropertiesManager::property_from_string( std::string const & property )
{
	return get_instance()->string_to_property_map().find( property )->second;
}

std::string const &
AtomPropertiesManager::string_from_property( AtomProperty const property )
{
	return get_instance()->property_to_string_map().find( property )->second;
}

// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
AtomPropertiesManager::AtomPropertiesManager() {}

}  // namespace chemical
}  // namespace core
