// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/AtomProperties.cc
/// @brief   Method definitions for AtomProperties.
/// @author  Labonte <JWLabonte@jhu.edu>

// Unit header
#include <core/chemical/AtomProperties.hh>
#include <core/chemical/AtomPropertiesManager.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>
#include <iostream>


// Construct tracer.
static THREAD_LOCAL basic::Tracer TR( "core.chemical.AtomProperties" );


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
AtomProperties::AtomProperties() : utility::pointer::ReferenceCount()
{
	init();
}

// Copy constructor
AtomProperties::AtomProperties( AtomProperties const & object_to_copy ) :
	utility::pointer::ReferenceCount( object_to_copy )
{
	copy_data( *this, object_to_copy );
}

// Destructor
AtomProperties::~AtomProperties() = default;

// Assignment operator
AtomProperties &
AtomProperties::operator=( AtomProperties const & object_to_copy )
{
	if ( this != &object_to_copy ) {
		copy_data( *this, object_to_copy );
	}
	return *this;
}

// Comparison operator
bool
AtomProperties::operator==( AtomProperties const & properties ) const
{
	return atom_property_status_ == properties.atom_property_status_;
}


// Standard Rosetta methods ///////////////////////////////////////////////////
// General methods
void
AtomProperties::show( std::ostream & output ) const
{
	using namespace std;

	output << "   Properties:";
	utility::vector1< string > const & properties( get_list_of_properties() );
	Size const n_properties( properties.size() );
	for ( uint i( 1 ); i <= n_properties; ++i ) {
		output << ' ' << properties[ i ];
	}
	output << endl;
}


// Accessors/Mutators
bool
AtomProperties::has_property( std::string const & property ) const
{
	return atom_property_status_[ AtomPropertiesManager::get_instance()->property_from_string( property ) ];
}

void
AtomProperties::set_property( std::string const & property, bool const setting )
{
	atom_property_status_[ AtomPropertiesManager::get_instance()->property_from_string( property ) ] = setting;
}


// Generate and return a list of strings representing the properties of this Atom.
utility::vector1< std::string >
AtomProperties::get_list_of_properties() const
{
	utility::vector1< std::string > list;

	for ( AtomProperty property = FIRST_ATOM_PROPERTY; property <= N_ATOM_PROPERTIES; ++property ) {
		if ( atom_property_status_[ property ] ) {
			list.push_back( AtomPropertiesManager::get_instance()->string_from_property( property ) );
		}
	}
	return list;
}


// Private methods ////////////////////////////////////////////////////////////
// Initialize data members.
void
AtomProperties::init()
{
	atom_property_status_.resize( N_ATOM_PROPERTIES, false );
}

// Copy all data members from <from> to <to>.
void
AtomProperties::copy_data( AtomProperties & to, AtomProperties const & from )
{
	to.atom_property_status_ = from.atom_property_status_;
}


// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that AtomProperties can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, AtomProperties const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

// This allows one to use a for loop with AtomProperty enum values.
// This is safe, because the AtomProperty values are set automatically before compiling.
AtomProperty &
operator++( AtomProperty & property )
{
	property = static_cast< AtomProperty >( static_cast< int >( property ) + 1 );
	return property;
}

}  // namespace chemical
}  // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AtomProperties::save( Archive & arc ) const {
	arc( CEREAL_NVP( atom_property_status_ ) ); // utility::vector1<_Bool>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AtomProperties::load( Archive & arc ) {
	arc( atom_property_status_ ); // utility::vector1<_Bool>
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AtomProperties );
CEREAL_REGISTER_TYPE( core::chemical::AtomProperties )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_AtomProperties )
#endif // SERIALIZATION
