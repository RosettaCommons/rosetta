// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/AtomProperties.hh
/// @brief   Declarations and simple accessor/mutator definitions for AtomProperties.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_AtomProperties_HH
#define INCLUDED_core_chemical_AtomProperties_HH

// Unit headers
#include <core/chemical/AtomProperties.fwd.hh>
#include <core/chemical/AtomPropertiesManager.fwd.hh>
#include <core/chemical/AtomProperty.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>

#ifdef WIN32
#include <core/types.hh>
#endif

namespace core {
namespace chemical {

/// @details This is a container class for storing properties specific to a ResidueType's atoms.
/// These properties belong to a particular ResidueType's Atoms; they do not belong to an AtomType.
/// chemical::Atoms store both AtomTypes and AtomProperties.
class AtomProperties : public utility::pointer::ReferenceCount {
public:  // Standard methods //////////////////////////////////////////////////
	/// @brief  Default constructor
	AtomProperties();

	/// @brief  Copy constructor
	AtomProperties( AtomProperties const & object_to_copy );

	// Destructor
	virtual ~AtomProperties();

	// Assignment operator
	AtomProperties & operator=( AtomProperties const & object_to_copy );

	// Comparison operator
	bool operator==( AtomProperties const & properties ) const;


public:  // Standard Rosetta methods //////////////////////////////////////////
	/// @brief  Generate string representation of AtomProperties for debugging purposes.
	virtual void show( std::ostream & output=std::cout ) const;


public:  // Accessors/Mutators ////////////////////////////////////////////////
	/// @brief  Get whether or not this Atom has the requested property.
	inline
	bool
	has_property( AtomProperty const property ) const
	{
		return atom_property_status_[ property ];
	}

	/// @brief  Get whether or not this Atom has the requested property by string.
	bool has_property( std::string const & property ) const;


	/// @brief  Set the status of the given property for this Atom.
	void
	set_property( AtomProperty const property, bool const setting )
	{
		atom_property_status_[ property ] = setting;
	}

	/// @brief  Set the status of the given property for this Atom by string.
	void set_property( std::string const & property, bool const setting );


public:  // Other public methods //////////////////////////////////////////////
	/// @brief  Generate and return a list of strings representing the properties of this Atom.
	utility::vector1< std::string > get_list_of_properties() const;


private:  // Private methods //////////////////////////////////////////////////
	// Initialize data members.
	void init();

	// Copy all data members from <from> to <to>.
	void copy_data( AtomProperties & to, AtomProperties const & from );


private:  // Private data /////////////////////////////////////////////////////
	// Storage of general atom properties.
	utility::vector1< bool > atom_property_status_;  // indexed by AtomProperty
};  // class AtomProperties


// Insertion operator (overloaded so that AtomProperties can be "printed" in PyRosetta).
std::ostream & operator<<( std::ostream & output, AtomProperties const & object_to_output );

// This allows one to use a for loop with AtomProperty enum values.
AtomProperty & operator++( AtomProperty & property );

}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_AtomProperties_HH
