// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_generator/ConstraintGeneratorFactory.cc
/// @brief  Implementation of the class for instantiating arbitrary ConstraintGenerators
///         from a string --> ConstraintGeneratorCreator map
/// @author Tom Linsky ( tlinsky at uw dot edu )

// Unit headers
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>

// Package headers
#include <protocols/constraint_generator/ConstraintGenerator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorCreator.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
namespace protocols {
namespace constraint_generator {

ConstraintGeneratorFactory::ConstraintGeneratorFactory():
	utility::SingletonBase< ConstraintGeneratorFactory >(),
	creator_map_()
{}

void
ConstraintGeneratorFactory::factory_register( ConstraintGeneratorCreatorOP creator )
{
	if ( creator_map_.find( creator->keyname() ) != creator_map_.end() ) {
		std::string const err_msg = "Factory Name Conflict: Two or more ConstraintGeneratorCreators registered with the name " + creator->keyname();
		utility_exit_with_message(  err_msg );
	}
	creator_map_[ creator->keyname() ] = creator;
}

bool ConstraintGeneratorFactory::has_type( std::string const & selector_type ) const
{
	return creator_map_.find( selector_type ) != creator_map_.end();
}

ConstraintGeneratorOP
ConstraintGeneratorFactory::new_constraint_generator(
	std::string const & constraint_generator_name,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) const
{
	if ( ! has_type( constraint_generator_name ) ) {
		std::string err_msg =  "No ConstraintGeneratorCreator with the name '" + constraint_generator_name + "' has been registered with the ConstraintGeneratorFactory";
		throw utility::excn::EXCN_Msg_Exception( err_msg );
	}
	auto iter = creator_map_.find( constraint_generator_name );
	ConstraintGeneratorOP new_constraint_generator = iter->second->create_constraint_generator();
	new_constraint_generator->parse_my_tag( tag, datamap );
	return new_constraint_generator;
}



void
ConstraintGeneratorFactory::define_constraint_generator_xml_schema_group( utility::tag::XMLSchemaDefinition & xsd ) const{
	try{
		utility::tag::define_xml_schema_group(
			creator_map_,
			constraint_generator_xml_schema_group_name(),
			& complex_type_name_for_constraint_generator,
			xsd );
	} catch( utility::excn::EXCN_Msg_Exception const & e ) {
		throw utility::excn::EXCN_Msg_Exception( "Could not generate an XML Schema for Constraints from ConstraintFactory; offending class"
			" must call protocols::constraint_generator::complex_type_name_for_constraint when defining"
			" its XML Schema\n" + e.msg() );
	}
}

std::string
ConstraintGeneratorFactory::constraint_generator_xml_schema_group_name(){
	return "constraint_generator";
}

std::string
ConstraintGeneratorFactory::complex_type_name_for_constraint_generator( std::string const & constraint_name ){
	return "constraint_generator_" + constraint_name + "_complex_type";
}

void
ConstraintGeneratorFactory::xsd_constraint_generator_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & constraint_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes)
{
	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_constraint_generator )
		.element_name( constraint_type )
		.description( description )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );
}

} //namespace constraint_generator
} //namespace protocols
