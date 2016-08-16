// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/AtomPropertiesManager.hh
/// @brief   Method declarations for AtomPropertiesManager.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_AtomPropertiesManager_HH
#define INCLUDED_core_chemical_AtomPropertiesManager_HH

// Unit headers
#include <core/chemical/AtomPropertiesManager.fwd.hh>
#include <core/chemical/AtomProperty.hh>

// Utility headers
#include <utility/SingletonBase.hh>

// C++ header
#include <string>
#include <map>


namespace core {
namespace chemical {

/// @details  This class is a singleton and manages AtomProperties enum mappings.
class AtomPropertiesManager: public utility::SingletonBase< AtomPropertiesManager > {
public:  // Declare friends ///////////////////////////////////////////////////
	friend class utility::SingletonBase< AtomPropertiesManager >;


public:  // Static constant data access ///////////////////////////////////////
	static AtomProperty const & property_from_string( std::string const & property );

	static std::string const & string_from_property( AtomProperty const property );


private:  // Private methods //////////////////////////////////////////////////
	// Empty constructor
	AtomPropertiesManager();

	// Singleton-creation function for use with utility::thread::threadsafe_singleton
	static AtomPropertiesManager * create_singleton_instance();

	// This private method is defined in AtomProperty_mappings.cc,
	// which will eventually be auto-generated.
	std::map< AtomProperty, std::string > const & property_to_string_map();

	// This private method is defined in AtomProperty_mappings.cc,
	// which will eventually be auto-generated.
	std::map< std::string, AtomProperty > const & string_to_property_map();


private:  // Private data /////////////////////////////////////////////////////
	std::map< AtomProperty, std::string > property_to_string_map_;
	std::map< std::string, AtomProperty > string_to_property_map_;
};

}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_AtomPropertiesManager_HH
