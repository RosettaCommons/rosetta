// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/RosettaScriptsSchemaValidator.cc
/// @brief  A singleton class for generating the schema accepted by the
///         RosettaScriptsParser and for holding the libxml2-library
///         schema validation object.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/rosetta_scripts/RosettaScriptsSchemaValidator.hh>

// Package Headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/parser/DataLoaderFactory.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <core/pack/palette/PackerPaletteFactory.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Basic headers
#include <basic/Tracer.hh>

// c++ headers
#include <map>
#include <set>
#include <sstream>

namespace protocols {
namespace rosetta_scripts {

static basic::Tracer TR( "protocols.rosetta_scripts.RosettaScriptsSchemaValidator" );

RosettaScriptsSchemaValidator::~RosettaScriptsSchemaValidator() {}

RosettaScriptsSchemaValidator::RosettaScriptsSchemaValidator() :
	validator_( new utility::tag::XMLValidator )
{
	using namespace utility::tag;

	std::ostringstream oss;
	oss << "If you are seeing this message, the internally-generated XML Schema"
		" for rosetta_scripts could not be properly generated\nThis failure occurred"
		" before your XML script that provided was examined. The error has been"
		" compiled into Rosetta and will need to be fixed by a developer.";
	std::string schema;
	try {
		schema = xsd_for_rosetta_scripts();
	} catch ( utility::excn::Exception const & e ) {
		oss << "An error was encountered while the string of the schema was being generated; this is before the schema is analyzed for whether it is correct or not.\n";
		oss << e.msg() << "\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  oss.str() );
	}

	XMLValidationOutput schema_valid_output;
	try {
		TR << "Initializing schema validator..." << std::endl;
		schema_valid_output = validator_->set_schema( schema );
		TR << "...done" << std::endl;
	} catch ( utility::excn::Exception const & e ) {
		oss << e.msg() << "\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  oss.str() );
	}
	if ( ! schema_valid_output.valid() ) {
		oss << "Global schema validation error: "
			" read the error message(s) below and fix your"
			" schema in the C++ code.\n\n";
		oss << "Errors: " << schema_valid_output.error_messages() << "\n";
		oss << "Warnings: " << schema_valid_output.warning_messages() << "\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  oss.str() );
	}
}

std::string
RosettaScriptsSchemaValidator::xsd_for_rosetta_scripts()
{
	TR << "Generating XML Schema for rosetta_scripts..." << std::endl;
	using namespace utility::tag;
	XMLSchemaDefinition xsd;

	write_ROSETTASCRIPTS_complex_type( xsd );

	// Finally, the root of the tree.
	// TO DO: separate the definition of the ROSETTA_SCRIPTS complex type and everything
	// that goes inside of it into its own function so that the ROSETTA_SCRIPTS element
	// becomes optional.
	XMLSchemaElement rosetta_scripts_element;
	rosetta_scripts_element.name( rosetta_scripts_element_name() )
		.type_name( rosetta_scripts_complex_type_naming_func( rosetta_scripts_element_name() ));
	xsd.add_top_level_element( rosetta_scripts_element );

	TR << "...done" << std::endl;

	return xsd.full_definition();

}

void
RosettaScriptsSchemaValidator::write_ROSETTASCRIPTS_complex_type( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	// early exit if we've been here already
	if ( xsd.has_top_level_element( rosetta_scripts_complex_type_naming_func( rosetta_scripts_element_name() ) ) ) return;

	// We have lots of complexTypes to define; we need them for the following elements:
	//
	// RESOURCES
	// MOVERS
	// FILTERS
	// All the data loaders, (TASKOPERATIONS, SCOREFXNS, etc.)
	// IMPORT
	// PROTOCOLS -- PROTOCOLS subelement is just the complexType for the ParsedProtocol mover!
	// OUTPUT
	// APPLY_TO_POSE ??
	// ROSETTASCRIPTS
	//
	// Then we need complexTypes defined from all of the DataLoaders -- we'll ask the
	// DataLoaderFactory for these.

	// Do the ROSETTASCRIPTS complex type first, so that if we recurse here, we can
	// quit (Note that the MoverFactory will recurse to the RosettaScriptsParser via
	// the ParsedProtocol mover)

	XMLSchemaSimpleSubelementList rosetta_scripts_resources_subelement;
	rosetta_scripts_resources_subelement.add_already_defined_subelement(
		"RESOURCES", & rosetta_scripts_complex_type_naming_func );

	XMLSchemaSimpleSubelementList rosetta_scripts_initial_subelements;
	for ( auto const & data_loader_pair : protocols::parser::DataLoaderFactory::get_instance()->loader_map() ) {
		rosetta_scripts_initial_subelements.add_already_defined_subelement(
			data_loader_pair.first, data_loader_pair.second->schema_ct_naming_function() );
	}
	rosetta_scripts_initial_subelements.add_already_defined_subelement(
		"MOVERS", & protocols::parser::DataLoaderFactory::data_loader_ct_namer );
	rosetta_scripts_initial_subelements.add_already_defined_subelement(
		"FILTERS", & protocols::parser::DataLoaderFactory::data_loader_ct_namer );
	rosetta_scripts_initial_subelements.add_already_defined_subelement(
		"APPLY_TO_POSE", & rosetta_scripts_complex_type_naming_func );
	rosetta_scripts_initial_subelements.add_already_defined_subelement(
		"IMPORT", & rosetta_scripts_complex_type_naming_func );

	XMLSchemaSimpleSubelementList rosetta_scripts_protocols_subelement;
	rosetta_scripts_protocols_subelement.add_already_defined_subelement_w_alt_element_name(
		"PROTOCOLS", ParsedProtocol::mover_name(), & moves::complex_type_name_for_mover );

	XMLSchemaSimpleSubelementList rosetta_scripts_output_subelement;
	rosetta_scripts_output_subelement.add_already_defined_subelement(
		"OUTPUT", & rosetta_scripts_complex_type_naming_func );

	XMLSchemaComplexTypeGenerator rosetta_scripts_ct;
	rosetta_scripts_ct.element_name( rosetta_scripts_element_name() )
		.complex_type_naming_func( & rosetta_scripts_complex_type_naming_func )
		.description( "The main ROSETTASCRIPTS block allows the definition of many different types of objects which can be generated"
		" from a text file. The combination of Movers and Filters can define an end-to-end protocol, and the use of the MultiplePoseMover"
		" can define pseudo-parallel diversification and pruning of a set of structures." )
		.add_ordered_subelement_set_as_optional( rosetta_scripts_resources_subelement )
		.add_ordered_subelement_set_as_repeatable( rosetta_scripts_initial_subelements )
		.add_ordered_subelement_set_as_required( rosetta_scripts_protocols_subelement )
		.add_ordered_subelement_set_as_optional( rosetta_scripts_output_subelement )
		.write_complex_type_to_schema( xsd );


	// Make sure the schemas for the DataLoaders, for all of the Movers, and all of the Filters are
	// written to the XSD.
	for ( auto const & data_loader_pair : protocols::parser::DataLoaderFactory::get_instance()->loader_map() ) {
		data_loader_pair.second->provide_xml_schema( xsd );
	}
	moves::MoverFactory::get_instance()->define_mover_xml_schema( xsd );
	filters::FilterFactory::get_instance()->define_filter_xml_schema( xsd );

	// RESOURCES
	AttributeList resource_attributes;
	resource_attributes + XMLSchemaAttribute::required_attribute( "name", xs_string, "The name of the resource that will"
		" be used in this XML file; to be used, this resource must have been declared in a resource-definition file and"
		" read in by the ResourceManager" );
	XMLSchemaSimpleSubelementList resources_subelements;
	resources_subelements.add_simple_subelement( "Resource", resource_attributes, "Each Resource to be used in a given job"
		" must be listed first at the beginning of the XML file so that the parser may request the Resource from the"
		" ResourceManager and put into the DataMap" );
	XMLSchemaComplexTypeGenerator resources_ctgen;
	resources_ctgen.element_name( "RESOURCES" )
		.complex_type_naming_func( & rosetta_scripts_complex_type_naming_func )
		.description( "Declare all of the Resources that should be loaded into the DataMap from the ResourceManager for this job; JD3 only." )
		.set_subelements_repeatable( resources_subelements )
		.write_complex_type_to_schema( xsd );

	// MOVERS
	XMLSchemaSimpleSubelementList movers_subelements;
	movers_subelements.add_group_subelement( & moves::MoverFactory::mover_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator movers_ct_gen;
	movers_ct_gen.element_name( "MOVERS" )
		.complex_type_naming_func( & protocols::parser::DataLoaderFactory::data_loader_ct_namer )
		.description( "The set of all of the Movers that are to be used in the script" )
		.set_subelements_repeatable( movers_subelements )
		.write_complex_type_to_schema( xsd );

	// FILTERS
	XMLSchemaSimpleSubelementList filters_subelements;
	filters_subelements.add_group_subelement( & filters::FilterFactory::filter_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator filters_ct_gen;
	filters_ct_gen.element_name( "FILTERS" )
		.complex_type_naming_func( & protocols::parser::DataLoaderFactory::data_loader_ct_namer )
		.description( "The set of all of the Filters that are to be used in the script" )
		.set_subelements_repeatable( filters_subelements )
		.write_complex_type_to_schema( xsd );

	// IMPORT
	AttributeList import_attributes;
	import_attributes
		+ XMLSchemaAttribute( "taskoperations", xsct_task_operation_comma_separated_list, "A comma separated list of TaskOperations that have been"
		" defined at a higher level than the MultiplePoseMover that this IMPORT statement is inside of" )
		+ XMLSchemaAttribute( "movers", xs_string, "A comma separated list of Movers that have been"
		" defined at a higher level than the MultiplePoseMover that this IMPORT statement is inside of" )
		+ XMLSchemaAttribute( "filters", xs_string, "A comma separated list of Filters that have been"
		" defined at a higher level than the MultiplePoseMover that this IMPORT statement is inside of" );
	XMLSchemaComplexTypeGenerator import_ct_gen;
	import_ct_gen.element_name( "IMPORT" )
		.complex_type_naming_func( & rosetta_scripts_complex_type_naming_func )
		.description( "The IMPORT statement is meant to be used in conjuction with the MultiplePoseMover. It allows"
		" users who have defined Movers, Filters, and TaskOperations at a higher level of the protocol to avoid"
		" having to redefine them in order to use them within the MultiplePoseMover. The classes that are imported"
		" will be re-parsed within the inner context. This parsing occurs before any of sub-elements of the"
		" ROSETTASCRIPTS block are loaded, reguardless of the order in which the IMPORT statement appears. Multiple"
		" IMPORT statements may be needed, e.g. if certain TaskOperations must be loaded before other Movers which"
		" require those TaskOperations are parsed." )
		.add_attributes( import_attributes )
		.write_complex_type_to_schema( xsd );

	// OUTPUT
	// Output element only contains attributes defined by parse_score_function
	XMLSchemaComplexTypeGenerator output_ct;
	AttributeList attributes;
	rosetta_scripts::attributes_for_parse_score_function( attributes );
	output_ct
		.complex_type_naming_func( & rosetta_scripts_complex_type_naming_func )
		.element_name( "OUTPUT" )
		.add_attributes( attributes )
		.description( "The OUTPUT element controls the final score function that is used by RosettaScripts before outputting a Pose."
		" If for example you have been using a score function with constraints in your protocol, but at the end, want to output the"
		" unconstrained scores, then you can set the score function in the OUTPUT element to some ScoreFunction that has constraint"
		" weights of zero and then RosettaScripts will rescore the Pose with that ScoreFunction right before writing the scores to"
		" disk." )
		.write_complex_type_to_schema( xsd );

	// APPLY_TO_POSE
	XMLSchemaSimpleSubelementList apply_to_pose_subelements;
	apply_to_pose_subelements.add_group_subelement( & moves::MoverFactory::mover_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator apply_to_pose_ct;
	apply_to_pose_ct.element_name( "APPLY_TO_POSE" )
		.complex_type_naming_func( & rosetta_scripts_complex_type_naming_func )
		.description( "The APPLY_TO_POSE block is deprecated and should not be used."
		" Any non-empty APPLY_TO_POSE block will now raise an error."
		" Move any movers from the APPLY_TO_POSE section into the PROTOCOLS block."
		" You may need to adjust some movers and filters to use reference poses and add a "
		" SavePoseMover to the PROTOCOLS block, right after the moved movers."
		)
		.set_subelements_repeatable( apply_to_pose_subelements )
		.write_complex_type_to_schema( xsd );

}


std::string
RosettaScriptsSchemaValidator::rosetta_scripts_element_name()
{
	return "ROSETTASCRIPTS";
}

std::string
RosettaScriptsSchemaValidator::rosetta_scripts_complex_type_naming_func( std::string const & element_name )
{
	return "rosetta_scripts_parser_" + element_name + "_type";
}


utility::tag::XMLSchemaValidatorCOP
RosettaScriptsSchemaValidator::validator() const
{
	return validator_;
}


} //namespace rosetta_scripts
} //namespace protocols


