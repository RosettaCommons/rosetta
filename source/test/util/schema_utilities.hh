// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/util/schema_utilities.hh
/// @brief  common utility functions for testing the quality of the documentation for XML schema
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_util_schema_utilities_HH
#define INCLUDED_util_schema_utilities_HH

#include <cxxtest/TestSuite.h>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/XMLSchemaValidation.hh>
#include <boost/function.hpp>
#include <map>
#include <string>

#include <basic/Tracer.hh>

static basic::Tracer TR2( "schema_utilities" );

inline
void recurse_through_subtags_for_attribute_descriptions(
	utility::tag::TagCOP tag,
	std::string const & residue_selector_name
) {
	for ( auto iter : tag->getTags() ) {
		if ( iter->getName() == "xs:attribute" ) {

			if ( ! iter->hasTag( "xs:annotation" ) ) {
				std::cerr << "Failed to find documentation in an \"xs:annotation\" sub-element for " << iter->getOption< std::string >( "name", "(unnamed attribute)" ) << " from complexType " << tag->getOption< std::string >( "name", "(unnamed complex type? this should never happen)" ) << std::endl;
				TS_ASSERT( iter->hasTag( "xs:annotation" ) );
				continue;
			}
			if ( ! iter->getTag( "xs:annotation" )->hasTag( "xs:documentation" ) ) {
				std::cerr << "Failed to find \"xs:documentation\" sub-element of an \"xs:annotation\" sub-element for " << iter->getOption< std::string >( "name", "(unnamed attribute)" ) << " from complexType " << tag->getOption< std::string >( "name", "(unnamed complex type? this should never happen)" ) << std::endl;
				TS_ASSERT( iter->getTag( "xs:annotation" )->hasTag( "xs:documentation" ) );
				continue;
			}
		} else {
			recurse_through_subtags_for_attribute_descriptions( iter, residue_selector_name );
		}
	}
}

// CreatorOP should define -> and then provide access to a class that defines
// provide_xml_schema taking a utility::tag::XMLSchemaDefinition &.
template < class CreatorOP >
void
ensure_all_cts_for_creators_have_documentation_strings(
	std::map< std::string, CreatorOP > const & creator_map,
	std::string const & name_of_base_class,
	boost::function< std::string ( std::string const & ) > ct_naming_function
)
{
	using namespace utility::tag;
	for ( auto iter : creator_map ) {
		XMLSchemaDefinition xsd;
		iter.second->provide_xml_schema( xsd );
		std::string full_def = xsd.full_definition();
		//TR2 << "full def: " << iter.first << "\n" << full_def << std::endl;
		TagCOP tag( Tag::create( full_def ) );

		// now we need to search for the xs:complexType with name = ct_naming_function( iter )
		std::string name = ct_naming_function( iter.first );
		TagCOP ct_tag;
		for ( auto child : tag->getTags() ) {
			if ( child->getName() == "xs:complexType" ) {
				ct_tag = child;
				break;
			}
		}
		if ( ! ct_tag ) {
			std::cerr << "Failed to identify the complex type in schema for " << name_of_base_class << " named " << iter.first << std::endl;
			TS_ASSERT( ct_tag );
			continue;
		}
		if ( ! ct_tag->hasTag( "xs:annotation" ) ) {
			std::cerr << "Failed to find class-level documentation in the form of an \"xs:annotation\" sub-element for " << name_of_base_class << " named " << iter.first << std::endl;
			TS_ASSERT( ct_tag->hasTag( "xs:annotation" ) );
			continue;
		}
		if ( ! ct_tag->getTag( "xs:annotation" )->hasTag( "xs:documentation" ) ) {
			std::cerr << "Failed to find class-level documentation in the form of an \"xs:documentation\" sub-element of an \"xs:annotation\" sub-element for " << name_of_base_class << " named " << iter.first << std::endl;
			TS_ASSERT( ct_tag->getTag( "xs:annotation" )->hasTag( "xs:documentation" ) );
			continue;
		}
	}
}


inline
void
ensure_rosetta_scripts_like_XSD_validates_w_group(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & group_name
)
{
	using namespace utility::tag;

	XMLSchemaComplexTypeOP rosetta_scripts_type( new XMLSchemaComplexType );
	XMLSchemaModelGroupOP rosetta_scripts_seq( new XMLSchemaModelGroup( xsmgt_sequence ));
	XMLSchemaElementOP rosetta_scripts_element( new XMLSchemaElement );

	XMLSchemaElementOP dummy_element( new XMLSchemaElement );
	//XMLSchemaModelGroupOP dummy_choice( new XMLSchemaModelGroup( xsmgt_choice ));
	XMLSchemaComplexTypeOP dummy_type( new XMLSchemaComplexType );

	XMLSchemaModelGroupOP group_ref( new XMLSchemaModelGroup( xsmgt_group ) );
	group_ref->group_name( group_name );
	group_ref->min_occurs( 0 );
	group_ref->max_occurs( xsminmax_unbounded );
	dummy_type->name( "DUMMY_Type" );
	dummy_type->set_model_group( group_ref );

	xsd.add_top_level_element( *dummy_type );

	dummy_element->name( "DUMMY" );
	dummy_element->type_name( "DUMMY_Type" );
	dummy_element->min_occurs( 0 );

	rosetta_scripts_seq->append_particle( dummy_element );
	rosetta_scripts_type->set_model_group( rosetta_scripts_seq );

	rosetta_scripts_element->name( "ROSETTASCRIPTS" );
	rosetta_scripts_element->element_type_def( rosetta_scripts_type );

	xsd.add_top_level_element( *rosetta_scripts_element );

	XMLValidationOutput output = test_if_schema_is_valid( xsd.full_definition() );
	TS_ASSERT( output.valid() );
	if ( ! output.valid() ) {
		TR2 << output.error_messages() << std::endl;
	}
}

inline
void
ensure_rosetta_scripts_like_XSD_validates_w_ct(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & ct_name
)
{
	using namespace utility::tag;

	XMLSchemaComplexTypeOP rosetta_scripts_type( new XMLSchemaComplexType );
	XMLSchemaModelGroupOP rosetta_scripts_seq( new XMLSchemaModelGroup( xsmgt_sequence ));
	XMLSchemaElementOP rosetta_scripts_element( new XMLSchemaElement );

	XMLSchemaElementOP dummy_element( new XMLSchemaElement );
	dummy_element->name( "DUMMY" );
	dummy_element->type_name( ct_name );
	dummy_element->min_occurs( 0 );

	rosetta_scripts_seq->append_particle( dummy_element );
	rosetta_scripts_type->set_model_group( rosetta_scripts_seq );

	rosetta_scripts_element->name( "ROSETTASCRIPTS" );
	rosetta_scripts_element->element_type_def( rosetta_scripts_type );

	xsd.add_top_level_element( *rosetta_scripts_element );

	XMLValidationOutput output = test_if_schema_is_valid( xsd.full_definition() );
	TS_ASSERT( output.valid() );
	if ( ! output.valid() ) {
		TR2 << output.error_messages() << std::endl;
	}
}

#endif
