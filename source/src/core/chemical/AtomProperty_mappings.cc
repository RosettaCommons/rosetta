// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/ResidueProperty_mappings.cc
/// @brief   Method definitions for private singleton methods declared in AtomPropertiesManager.hh.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @note    TODO: Auto-generate this file in the same way ResiduePropertys are generated.

// Unit header
#include <core/chemical/AtomPropertiesManager.hh>

// C++ header
#include <map>


namespace core {
namespace chemical {

// This private method is declared in AtomPropertiesManager.hh.
// It will eventually be auto-generated.
std::map< AtomProperty, std::string > const &
AtomPropertiesManager::property_to_string_map()
{
	using namespace std;

	if ( property_to_string_map_.empty() ) {
		property_to_string_map_.insert( make_pair( NO_ATOM_PROPERTY, "NO_ATOM_PROPERTY" ) );
		property_to_string_map_.insert( make_pair( H_DONOR, "H_DONOR" ) );
		property_to_string_map_.insert( make_pair( H_ACCEPTOR, "H_ACCEPTOR" ) );
		property_to_string_map_.insert( make_pair( POLAR_HYDROGEN, "POLAR_HYDROGEN" ) );
		property_to_string_map_.insert( make_pair( AROMATIC_HYDROGEN, "AROMATIC_HYDROGEN" ) );
		property_to_string_map_.insert( make_pair( HAS_ORBITALS, "HAS_ORBITALS" ) );
		property_to_string_map_.insert( make_pair( VIRTUAL_ATOM, "VIRTUAL_ATOM" ) );
		property_to_string_map_.insert( make_pair( REPULSIVE, "REPULSIVE" ) );
		property_to_string_map_.insert( make_pair( AROMATIC_CARBON_WITH_FREE_VALENCE, "AROMATIC_CARBON_WITH_FREE_VALENCE" ) );
	}
	return property_to_string_map_;
}

// This private method is declared in AtomPropertiesManager.hh.
// It will eventually be auto-generated.
std::map< std::string, AtomProperty > const &
AtomPropertiesManager::string_to_property_map()
{
	using namespace std;

	if ( string_to_property_map_.empty() ) {
		string_to_property_map_.insert( make_pair( "NO_ATOM_PROPERTY", NO_ATOM_PROPERTY ) );
		string_to_property_map_.insert( make_pair( "H_DONOR", H_DONOR ) );
		string_to_property_map_.insert( make_pair( "H_ACCEPTOR", H_ACCEPTOR ) );
		string_to_property_map_.insert( make_pair( "POLAR_HYDROGEN", POLAR_HYDROGEN ) );
		string_to_property_map_.insert( make_pair( "AROMATIC_HYDROGEN", AROMATIC_HYDROGEN ) );
		string_to_property_map_.insert( make_pair( "HAS_ORBITALS", HAS_ORBITALS ) );
		string_to_property_map_.insert( make_pair( "VIRTUAL_ATOM", VIRTUAL_ATOM ) );
		string_to_property_map_.insert( make_pair( "REPULSIVE", REPULSIVE ) );
		string_to_property_map_.insert( make_pair( "AROMATIC_CARBON_WITH_FREE_VALENCE", AROMATIC_CARBON_WITH_FREE_VALENCE ) );
	}
	return string_to_property_map_;
}

}  // namespace chemical
}  // namespace core
