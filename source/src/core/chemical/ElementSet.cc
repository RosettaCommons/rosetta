// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/ElementSet.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <core/chemical/ElementSet.hh>
#include <core/chemical/Element.hh>

// Project headers
#include <basic/Tracer.hh>
#include <utility>
#include <utility/io/izstream.hh>

// C++ headers
#include <fstream>
#include <string>

#ifdef WIN32
#include <utility/string_util.hh>
#endif

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {


static basic::Tracer tr( "core.chemical.ElementSet" );

ElementSet::ElementSet(std::string const & name): name_(name) {}

ElementSet::~ElementSet() = default;

/// @details Initialize an ElementSet from an external file "filename",
/// and set parameters and properties for each Element.
/// Refer to rosetta_database/chemical/elements/element_properties.txt
/// for file format
void
ElementSet::read_file( std::string const & filename )
{
	utility::io::izstream data( filename.c_str() );

	if ( !data.good() ) utility_exit_with_message( "Unable to open element file: "+filename );

	// now parse the rest of the file
	{
		using namespace basic;

		std::string line;
		while ( data.good() ) {
			char next;
			while ( static_cast<void>( next = data.peek() ) /*This static cast is necessary to silence a warning in the xcode compilation*/, next == ' ' || next == '\n' || next == '\t' ) { data.get(); } // Discard leading whitespace
			if ( data.peek() == '#' ) {
				getline( data,line ); // Discard the comment line
				continue;
			}
			if ( data.peek() == 'E' ) { // Beginning of "Element:" tag
				ElementOP element( new Element );
				data >> *element;
				if ( data.good() ) {
					elements_.push_back( element );
					std::string symbol( element->get_chemical_symbol() );
					if ( element_index_.count( symbol ) ) {
						if ( symbol != "Z" && symbol != "X" ) {
							// Z and X are used as dummy element symbols -- ignore
							utility_exit_with_message("ElementSet:: duplicate element symbol "+symbol);
						}
					} else {
						element_index_[ symbol ] = elements_.size();
					}
					tr.Debug << "New element: " << symbol << std::endl;
				}
			} else {
				getline( data, line );
				if ( data.good() ) {
					utility_exit_with_message("Problem with Elements file. Expecting 'Element:' tag, found:'" + line + "'");
				}
			}
		}
	}
}

/// @brief Check if there is an element_type associated with an element_symbol string
bool
ElementSet::contains_element_type( std::string const & element_symbol ) const
{
	std::locale loc;
	auto
		iter( element_index_.find( element_symbol ) );
	// If we can't find it straight-away, we may need to title case it.
	// (Cl versus CL)
	if ( iter == element_index_.end() && element_symbol.size() >=1 && element_symbol.size() <= 2 ) {
		std::string title( 1, std::toupper(element_symbol[0], loc) );
		if ( element_symbol.size() == 2 ) { title += std::tolower(element_symbol[1], loc); }
		iter = element_index_.find( title );
	}
	return iter != element_index_.end();
}


/// @brief Lookup the element index by the element enum
Size
ElementSet::element_index( core::chemical::element::Elements ele ) const
{
	// Probably not the most efficient way to do things, but ...
	return element_index( core::chemical::element::name_from_elements( ele ) );
}


/// @brief Lookup the element index by the element_symbol string
Size
ElementSet::element_index( std::string const & element_symbol ) const
{
	std::locale loc;
	auto
		iter( element_index_.find( element_symbol ) );
	// If we can't find it straight-away, we may need to title case it.
	// (Cl versus CL)
	if ( iter == element_index_.end() && element_symbol.size() >=1 && element_symbol.size() <= 2 ) {
		std::string title( 1, std::toupper(element_symbol[0], loc) );
		if ( element_symbol.size() == 2 ) { title += std::tolower(element_symbol[1], loc); }
		iter = element_index_.find( title );
	}
	if ( iter == element_index_.end() ) {
		utility_exit_with_message( "unrecognized element_symbol '"+element_symbol+"'" );
	}
	return iter->second;
}

/// @brief Lookup the element index by the element enum
ElementCOP
ElementSet::element( core::chemical::element::Elements ele ) const
{
	return elements_[ element_index( ele ) ];
}

/// @brief Lookup the element index by the element_symbol string
ElementCOP
ElementSet::element( std::string const & element_symbol ) const
{
	return elements_[ element_index( element_symbol ) ];
}


/// @brief Lookup an Element by 1-based indexing
ElementCOP
ElementSet::operator[] ( Size const index ) const
{
	return elements_[ index ];
}


} // chemical
} // core

#ifdef SERIALIZATION
#include <core/chemical/ChemicalManager.hh>

template < class Archive >
void core::chemical::serialize_element_set( Archive & arc, core::chemical::ElementSetCOP ptr )
{
	if ( ! ptr ) {
		bool ptr_is_nonnull( false );
		arc( CEREAL_NVP( ptr_is_nonnull ) );
	} else {
		bool ptr_is_nonnull( true );
		arc( CEREAL_NVP( ptr_is_nonnull ) );
		std::string typeset_name( ptr->name() ); // Assumes that the name can be used to extract it from the ChemicalManager
		arc( CEREAL_NVP( typeset_name ) );
	}
}
INSTANTIATE_FOR_OUTPUT_ARCHIVES( void, core::chemical::serialize_element_set, core::chemical::ElementSetCOP );

template < class Archive >
void core::chemical::deserialize_element_set( Archive & arc, core::chemical::ElementSetCOP & ptr )
{
	bool ptr_is_nonnull( true ); arc( ptr_is_nonnull );
	if ( ptr_is_nonnull ) {
		std::string typeset_name;
		arc( typeset_name );
		ptr = core::chemical::ChemicalManager::get_instance()->element_set( typeset_name );
	} else {
		ptr = nullptr;
	}
}
INSTANTIATE_FOR_INPUT_ARCHIVES( void, core::chemical::deserialize_element_set, core::chemical::ElementSetCOP & );

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_ElementSet )
#endif // SERIALIZATION
