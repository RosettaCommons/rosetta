// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/utility/tag/xml_schema_group_initialization.hh
/// @brief  a templated function for the initialization of xml schema groups
///         designed to be used by factories that contain maps of pointers to
///         WidgetCreators.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_tag_xml_schema_group_initialization_hh
#define INCLUDED_utility_tag_xml_schema_group_initialization_hh

// Unit headers
#include <utility/tag/XMLSchemaGeneration.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/backtrace.hh> // for debug_assert
#include <utility/excn/Exceptions.hh>

// Boost headers
#include <boost/function.hpp>

// C++ headers
#include <map>
#include <string>

namespace utility {
namespace tag {

/// @details Creator must define two functions: keyname() and provide_xml_schema()
///
/// @throws In the event that the Creator does not define a complex type with the
/// name coming from the complex_type_name_for_widget_func, this function will throw
/// a utility::excn::Excn_Msg_Excption.  Classes calling this function should wrap
/// it in a try/catch block and append extra information informing the user (presumably
/// the programmer) which function should be called when defining the schema for the
/// offending class.
template < class Creator >
void define_xml_schema_group(
	typename std::map< std::string, utility::pointer::shared_ptr< Creator > > const & creator_map,
	std::string const & widget_group_name,
	boost::function< std::string( std::string const & ) > const & complex_type_name_for_widget_func,
	XMLSchemaDefinition & xsd
)
{
	typedef std::map< std::string, utility::pointer::shared_ptr< Creator > > CreatorMap;

	if ( xsd.has_top_level_element( widget_group_name ) ) {
		// early exit to prevent infinite recursion and the duplication of effort
		return;
	}

	XMLSchemaModelGroup widget_group;
	widget_group.group_name( widget_group_name );
	XMLSchemaModelGroupOP widget_group_choice( new XMLSchemaModelGroup );
	widget_group_choice->type( xsmgt_choice );
	for ( typename CreatorMap::const_iterator iter = creator_map.begin(), iter_end = creator_map.end(); iter != iter_end; ++iter ) {
		XMLSchemaElementOP widget_element( new XMLSchemaElement );
		widget_element->name( iter->first ).type_name( complex_type_name_for_widget_func( iter->first ));
		widget_group_choice->append_particle( widget_element );
	}
	widget_group.append_particle( widget_group_choice );

	xsd.add_top_level_element( widget_group );

	// Now iterate across all Creators in the map, and have each one write their definition to the XSD
	for ( typename CreatorMap::const_iterator iter = creator_map.begin(), iter_end = creator_map.end(); iter != iter_end; ++iter ) {

		try {
			iter->second->provide_xml_schema( xsd );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			throw utility::excn::EXCN_Msg_Exception( "define_xml_schema_group: failed to define an xml schema for \"" + iter->first + "\"; message was:\n" + e.msg() );
		}

		if ( ! xsd.has_top_level_element( complex_type_name_for_widget_func( iter->first ) ) ) {
			//std::cout << "schema:\n" << xsd.full_definition() << std::endl;
			throw utility::excn::EXCN_Msg_Exception( "define_xml_schema_group: failed to detect a complex type of name \"" +
				complex_type_name_for_widget_func( iter->first ) + "\" for \"" +
				iter->first + "\"\n" );
		}
	}

}

}
}

#endif

