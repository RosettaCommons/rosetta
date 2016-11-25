// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/rdf/RDFFunctionFactory.cc
/// @brief  Factory for creating RDFBase objects
/// @author Sam DeLuca

// Unit Headers
#include <protocols/ligand_docking/rdf/RDFFunctionFactory.hh>
#include <protocols/ligand_docking/rdf/RDFFunctionCreator.hh>
#include <protocols/ligand_docking/rdf/RDFBase.hh>

// Package Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector0.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
// C++ Headers
#include <sstream>

//Auto Headers
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {
namespace rdf {

using std::endl;
using std::string;
using std::stringstream;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using protocols::filters::Filters_map;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;
using utility::tag::TagCOP;

static THREAD_LOCAL basic::Tracer tr( "protocols.ligand_docking.rdf.RDFFunctionFactory" );

/// @details Private constructor insures correctness of singleton.
RDFFunctionFactory::RDFFunctionFactory() {}

RDFFunctionFactory::~RDFFunctionFactory() {}

void
RDFFunctionFactory::factory_register(
	RDFFunctionCreatorCOP creator
) {
	types_[ creator->type_name() ] = creator;
}


RDFBaseOP
RDFFunctionFactory::get_rdf_function(std::string const & type_name)
{
	tr.Trace << "generate RDF function of type " << type_name << std::endl;
	RDFFunctionCreatorMap::const_iterator iter = types_.find( type_name );
	if ( iter != types_.end() ) {
		return iter->second->create_rdf_function();
	} else {
		stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized RDF Function "
			<< "'" << type_name << "'." << endl
			<< "check spelling or "
			<< "register a new RDF Function in the RDFFunctionFactory" << endl
			<< "known RDF Function types are:" << endl;

		for ( auto const & type : types_ ) {
			error_msg << "\t" << type.first << endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}
utility::vector1<std::string> RDFFunctionFactory::get_all_function_names()
{
	utility::vector1<std::string> collection;
	RDFFunctionCreatorMap::const_iterator iter = types_.begin(), end = types_.end();
	while ( iter != end ) {
		collection.push_back(iter->first);
		++iter;
	}
	return collection;

}
RDFBaseOP
RDFFunctionFactory::get_rdf_function(
	TagCOP const tag,
	basic::datacache::DataMap & data
)
{
	/*assert(tag->getName() == "RDF");

	string type_name;
	if ( !tag->hasOption("name") ) {
	utility_exit_with_message("'RDF' tags require a name field");
	} else {
	*/
	std::string type_name;
	type_name = tag->getName();

	RDFBaseOP rdf_function(get_rdf_function(type_name));

	rdf_function->parse_my_tag(tag, data);
	return rdf_function;
}

void
RDFFunctionFactory::define_rdf_function_group( utility::tag::XMLSchemaDefinition & xsd ){
	try {
		utility::tag::define_xml_schema_group(
			types_,
			rdf_function_group_name(),
			& rdf_function_ct_namer,
			xsd );
	} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
		throw utility::excn::EXCN_Msg_Exception( "Could not generate an XML Schema for RDFFunction from RDFFunctionFactory; offending class"
			" must call protocols::ligand_docking::rdf::rdf_function_ct_namer when defining"
			" its XML Schema\n" + e.msg() );
	}


}
std::string
RDFFunctionFactory::rdf_function_group_name(){
	return "rdf_function";
}
std::string
RDFFunctionFactory::rdf_function_ct_namer( std::string tag_name ){
	return "rdf_function_" + tag_name + "_complex_type";
}

void
RDFFunctionFactory::xsd_type_definition_w_attributes( utility::tag::XMLSchemaDefinition & xsd, std::string name, utility::tag::AttributeList & attlist, std::string description )
{
	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & rdf_function_ct_namer )
		.element_name( name )
		.add_optional_name_attribute( "Optional identifier for this RDFFunction" )
		.description( description )
		.add_attributes( attlist )
		.write_complex_type_to_schema( xsd );



}



} // namespace
} // namespace
} // namespace
