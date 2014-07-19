// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/ResidueProperties.cc
/// @brief   Method definitions for ResidueProperties.
/// @author  Labonte <JWLabonte@jhu.edu>

// Unit header
#include <core/chemical/ResidueProperties.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>
#include <iostream>


// Construct tracer.
static basic::Tracer TR( "core.chemical.ResidueProperties" );


namespace core {
namespace chemical {

using namespace core;

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
ResidueProperties::ResidueProperties() : utility::pointer::ReferenceCount()
{
	init();
}

// Copy constructor
ResidueProperties::ResidueProperties( ResidueProperties const & object_to_copy ) :
		utility::pointer::ReferenceCount( object_to_copy )
{
	copy_data( *this, object_to_copy );
}

// Assignment operator
ResidueProperties &
ResidueProperties::operator=( ResidueProperties const & object_to_copy )
{
	// Abort self-assignment.
	if ( this == &object_to_copy ) {
		return *this;
	}

	copy_data( *this, object_to_copy );
	return *this;
}

// Destructor
ResidueProperties::~ResidueProperties() {}


// Standard Rosetta methods ///////////////////////////////////////////////////
// General methods
void
ResidueProperties::show( std::ostream & output ) const
{
	using namespace std;
	using namespace utility;

	vector1< string > const & properties( get_list_of_properties() );
	Size const n_properties( properties.size() );

	output << " Properties:";
	for ( core::uint i = 1; i <= n_properties; ++i ) {
		output << ' ' << properties[ i ];
	}
	output << endl;
}


// Accessors/Mutators
bool
ResidueProperties::has_property( std::string const & property ) const
{
	// TODO: Wrap in try/catch.
	// tr.Warning << "WARNING:: unrecognized residue type property: " << property << std::endl;
	return general_property_status_[ get_property_from_string( property ) ];
}

void
ResidueProperties::set_property( std::string const & property, bool const setting )
{
	general_property_status_[ get_property_from_string( property ) ] = setting;
}


/// @note  This function was copied from its old location in ResidueType.
bool
ResidueProperties::has_variant_type( VariantType const & variant_type ) const
{
	using namespace std;

	return ( find( variant_types_.begin(), variant_types_.end(), variant_type ) != variant_types_.end() );
}

/// @note  This function was copied from its old location in ResidueType.
void
ResidueProperties::add_variant_type( VariantType const & variant_type )
{
	if ( ! has_variant_type( variant_type ) ) {
		variant_types_.push_back( variant_type );
	}
}


/// @note  This function was copied from its old location in ResidueType.
void
ResidueProperties::add_numeric_property( std::string const & tag, core::Real value )
{
	using namespace std;

	numeric_properties_.insert( make_pair( tag, value ) );
}

/// @note  This function was copied from its old location in ResidueType.
void
ResidueProperties::add_string_property( std::string const & tag, std::string value )
{
	using namespace std;

	string_properties_.insert( make_pair( tag, value ) );
}


// Generate and return a list of strings representing the properties of this ResidueType.
utility::vector1< std::string >
ResidueProperties::get_list_of_properties() const
{
	using namespace std;
	using namespace utility;

	vector1< string > list;

	for ( ResidueProperty property = FIRST_PROPERTY; property <= N_RESIDUE_PROPERTIES; ++property ) {
		if ( general_property_status_[ property ] ) {
			list.push_back( get_string_from_property( property ) );
		}
	}

	return list;
}


// Private methods ////////////////////////////////////////////////////////////
// Initialize data members.
void
ResidueProperties::init()
{
	general_property_status_.resize( N_RESIDUE_PROPERTIES, false );
}

// Copy all data members from <from> to <to>.
void
ResidueProperties::copy_data( ResidueProperties & to, ResidueProperties const & from )
{
	to.general_property_status_ = from.general_property_status_;
	to.variant_types_ = from.variant_types_;
	to.numeric_properties_ = from.numeric_properties_;
	to.string_properties_ = from.string_properties_;
}


// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that ResidueProperties can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, ResidueProperties const & object_to_output )
{
	object_to_output.show(output);
	return output;
}

// This allows one to use a for loop with ResidueProperty enum values.
// This is safe, because the ResidueProperty values are set automatically before compiling.
ResidueProperty &
operator++( ResidueProperty & property )
{
	property = static_cast< ResidueProperty >( static_cast< int >( property ) + 1 );
	return property;
}

}  // namespace chemical
}  // namespace core


