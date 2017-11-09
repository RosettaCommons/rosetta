// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   LegacyRequirementFactory.cc
/// @brief  Factory for creating Requirements objects
/// @author Tim Jacobs

// Unit Headers
#include <protocols/legacy_sewing/sampling/LegacyAssemblyMover.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementFactory.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyGlobalRequirement.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyIntraSegmentRequirement.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementCreator.hh>

// Package Headers
#include <basic/Tracer.hh>

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>

// C++ Headers
#include <sstream>

//Utility Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>



namespace protocols {
namespace legacy_sewing  {
namespace sampling {
namespace requirements {

static THREAD_LOCAL basic::Tracer tr("protocols.legacy_sewing.sampling.requirements.LegacyRequirementFactory");

/// @details Private constructor insures correctness of singleton.
LegacyRequirementFactory::LegacyRequirementFactory() {}

LegacyRequirementFactory::~LegacyRequirementFactory() {}

void
LegacyRequirementFactory::factory_register(
	LegacyGlobalRequirementCreatorCOP creator
) {
	global_types_[ creator->type_name() ] = creator;
}

void
LegacyRequirementFactory::factory_register(
	LegacyIntraSegmentRequirementCreatorCOP creator
) {
	intra_segment_types_[ creator->type_name() ] = creator;
}


LegacyGlobalRequirementOP
LegacyRequirementFactory::get_global_requirement(
	std::string const & type_name
) {
	tr.Trace << "Generating global requirement of type " << type_name << std::endl;
	LegacyGlobalRequirementCreatorMap::const_iterator iter = global_types_.find( type_name );
	if ( iter != global_types_.end() ) {
		return iter->second->create_requirement();
	} else {
		std::stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized LegacyGlobalRequirement "
			<< "'" << type_name << "'." << std::endl
			<< "check spelling or "
			<< "register a new Requirement in the LegacyRequirementFactory" << std::endl
			<< "known Requirement types are:" << std::endl;

		for ( auto const & type : global_types_ ) {
			error_msg << "\t" << type.first << std::endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}


LegacyIntraSegmentRequirementOP
LegacyRequirementFactory::get_intra_segment_requirement(
	std::string const & type_name
) {
	tr.Trace << "Generating intra-segment requirement of type " << type_name << std::endl;
	LegacyIntraSegmentRequirementCreatorMap::const_iterator iter = intra_segment_types_.find( type_name );
	if ( iter != intra_segment_types_.end() ) {
		return iter->second->create_requirement();
	} else {
		std::stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized LegacyIntraSegmentRequirement "
			<< "'" << type_name << "'." << std::endl
			<< "check spelling or "
			<< "register a new Requirement in the LegacyRequirementFactory" << std::endl
			<< "known Requirement types are:" << std::endl;

		for ( auto const & type : intra_segment_types_ ) {
			error_msg << "\t" << type.first << std::endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}


void
LegacyRequirementFactory::define_global_requirements_subelement( utility::tag::XMLSchemaDefinition & xsd ) const{
	using namespace utility::tag;

	define_xml_schema_group( global_types_, legacy_global_requirements_group_name(), & legacy_global_requirements_ct_namer, xsd );

	XMLSchemaSimpleSubelementList subs;
	subs.add_group_subelement( & legacy_global_requirements_group_name );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & LegacyAssemblyMover::legacy_assembly_mover_ct_namer )
		.element_name( "GlobalRequirements" )
		.description( "Subelement used to define global assembly requirements" )
		.set_subelements_repeatable( subs)
		.write_complex_type_to_schema( xsd );
}

void
LegacyRequirementFactory::define_intra_segment_requirements_subelement( utility::tag::XMLSchemaDefinition & xsd ) const{
	using namespace utility::tag;

	define_xml_schema_group( intra_segment_types_, legacy_intra_segment_requirements_group_name(), & legacy_intra_segment_requirements_ct_namer, xsd );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "index", xsct_non_negative_integer, "Which segment to apply these requirements to" );

	XMLSchemaSimpleSubelementList subs;
	subs.add_group_subelement( & legacy_intra_segment_requirements_group_name );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & LegacyAssemblyMover::legacy_assembly_mover_ct_namer )
		.element_name( "IntraSegmentRequirements" )
		.description( "Subelement used to define intra-segment assembly requirements" )
		.add_attributes( attlist )
		.set_subelements_repeatable( subs)
		.write_complex_type_to_schema( xsd );
}

std::string LegacyRequirementFactory::legacy_global_requirements_ct_namer( std::string const & input ){
	return "LegacyGlobalRequirements_" + input + "_type";
}
std::string LegacyRequirementFactory::legacy_global_requirements_group_name(){
	return "GlobalRequirements";
}
std::string LegacyRequirementFactory::legacy_intra_segment_requirements_ct_namer( std::string const & input ){
	return "LegacyIntraSegRequirements_" + input + "_type";
}
std::string LegacyRequirementFactory::legacy_intra_segment_requirements_group_name(){
	return "IntraSegmentRequirements";
}

//utility::vector1<std::string> LegacyRequirementFactory::get_all_features_names()
//{
// utility::vector1<std::string> collection;
// LegacyRequirementCreatorMap::const_iterator iter = types_.begin();
// while ( iter != types_.end() ) {
//  collection.push_back(iter->first);
//  iter++;
// }
// return collection;
//}

} // namespace
} // namespace
} // namespace
} // namespace
