// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/CompoundArchitect.cc
/// @brief Architect that creates a StructureData using multiple architects
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/architects/CompoundArchitect.hh>
#include <protocols/denovo_design/architects/CompoundArchitectCreator.hh>
// Protocol headers
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.hh>
#include <protocols/denovo_design/architects/PoseArchitect.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>
#include <protocols/denovo_design/connection/util.hh>

// Basic/Utililty headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.architects.CompoundArchitect" );

namespace protocols {
namespace denovo_design {
namespace architects {

CompoundArchitect::CompoundArchitect( std::string const & id_value ):
	DeNovoArchitect( id_value ),
	architects_(),
	connections_()
{
}

CompoundArchitect::~CompoundArchitect()
{}

CompoundArchitect::DeNovoArchitectOP
CompoundArchitect::clone() const
{
	return DeNovoArchitectOP( new CompoundArchitect( *this ) );
}

std::string
CompoundArchitect::type() const
{
	return CompoundArchitect::class_name();
}



void
CompoundArchitectCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	CompoundArchitect::provide_xml_schema( xsd );
}

std::string
CompoundArchitect::subelement_ct_namer( std::string tag_name ){
	return "compound_architect_subelement_" + tag_name + "_complex_type";
}
/*
std::string
CompoundArchitect::architect_subelement_ct_namer( std::string tag_name ){
return "compound_architect_architect_subtag_" + tag_name + "_complex_type";
}
std::string
CompoundArchitect::connection_subelement_ct_namer( std::string tag_name ){
return "compound_architect_connection_subtag_" + tag_name + "_complex_type";
}
*/
std::string
CompoundArchitect::pairing_subelement_ct_namer( std::string tag_name ) {
	return "pairing_architect_connection_subtag_" + tag_name + "_complex_type";
}

std::string
CompoundArchitect::pairing_group_namer(){
	return "compound_architect_pairing_group";
}
void
CompoundArchitect::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	//No attributes, just subtags
	using namespace utility::tag;
	//Possible subelements = "Architects", "Connections", "Pairing"
	DeNovoArchitectFactory::get_instance()->define_architect_group( xsd );
	//None of the subelements have any attributes, just more subelements
	XMLSchemaSimpleSubelementList architect_subelements;
	architect_subelements.add_group_subelement( & DeNovoArchitectFactory::architect_group_name );

	XMLSchemaComplexTypeGenerator architect_ct_gen;
	architect_ct_gen.element_name( "Architects" )
		.set_subelements_repeatable( architect_subelements )
		.complex_type_naming_func( & subelement_ct_namer )
		.description( "Use to build backbone from a combination of other architects" )
		.write_complex_type_to_schema( xsd );

	AttributeList connections_subelement_attributes;
	connections_subelement_attributes
		+ XMLSchemaAttribute::required_attribute( "name", xs_string, "XRW TO DO" );
	XMLSchemaSimpleSubelementList connections_subelements;
	connections_subelements
		.add_simple_subelement( "Connection", connections_subelement_attributes, "XRW TO DO" );
	XMLSchemaComplexTypeGenerator connections_ct_gen;
	connections_ct_gen.element_name( "Connections" )
		.set_subelements_repeatable( connections_subelements )
		.complex_type_naming_func( & subelement_ct_namer )
		.description( "XRW TO DO" )
		.write_complex_type_to_schema( xsd );

	define_pairing_group( xsd );
	XMLSchemaSimpleSubelementList pairing_subelements;
	pairing_subelements.add_group_subelement( & pairing_group_namer );
	//all of the SegmentPairingOP derived classes are valid
	XMLSchemaComplexTypeGenerator pairing_ct_gen;
	pairing_ct_gen.element_name( "Pairing" )
		.set_subelements_repeatable( pairing_subelements )
		.complex_type_naming_func( & subelement_ct_namer )
		.description( "XRW TO DO" )
		.write_complex_type_to_schema( xsd );



	AttributeList attlist;
	DeNovoArchitect::add_common_denovo_architect_attributes( attlist );

	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_already_defined_subelement( "Pairing",  & subelement_ct_namer )
		.add_already_defined_subelement( "Connections", & subelement_ct_namer )
		.add_already_defined_subelement( "Architects", & subelement_ct_namer );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( class_name() )
		.complex_type_naming_func( & DeNovoArchitectFactory::complex_type_name_for_architect )
		.set_subelements_repeatable( subelements )
		.description( "XRW TO DO" )
		.write_complex_type_to_schema( xsd );
}



void
CompoundArchitect::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	for ( utility::tag::Tag::tags_t::const_iterator t=tag->getTags().begin(); t!=tag->getTags().end(); ++t ) {
		if ( (*t)->getName() == "Architects" ) parse_architect_tags( *t, data );
		else if ( (*t)->getName() == "Connections" ) parse_connection_tags( *t, data );
		else if ( (*t)->getName() == "Pairing" ) parse_pairing_tags( **t );
		else {
			std::stringstream msg;
			msg << type() << "::parse_tag(): The provided subtag name (" << (*t)->getName()
				<< ") is not valid. Valid subtag names for " << type() << " are: "
				<< "Architects, Connections, Pairing." << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
		}
	}
}

void
CompoundArchitect::parse_architect_tags( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	for ( utility::tag::Tag::tags_t::const_iterator t=tag->getTags().begin(); t!=tag->getTags().end(); ++t ) {
		parse_architect_tag( *t, data );
	}
}

void
CompoundArchitect::parse_architect_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	DeNovoArchitectOP new_architect = DeNovoArchitectFactory::get_instance()->create_from_tag( tag, data );
	architects_.push_back( new_architect );
}

void
CompoundArchitect::parse_connection_tags( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	for ( utility::tag::Tag::tags_t::const_iterator t=tag->getTags().begin(); t!=tag->getTags().end(); ++t ) {
		parse_connection_tag( *t, data );
	}
}


void
CompoundArchitect::define_pairing_group( utility::tag::XMLSchemaDefinition & xsd ){
	//components::StrandPairing::class_name()
	//components::HelixPairing::class_name()
	//components::HelixSheetPairing::class_name()
	using namespace utility::tag;

	AttributeList strand_pairing_attlist;

	//Define the complex types for all of these
	XMLSchemaComplexTypeGenerator strand_pairing_ct;
	strand_pairing_ct.element_name( components::StrandPairing::class_name() )
		.complex_type_naming_func ( & pairing_subelement_ct_namer )
		.description( "XRW TO DO" )
		.add_attributes( strand_pairing_attlist )
		.write_complex_type_to_schema( xsd );


	//random and null have no attributes
	XMLSchemaComplexTypeGenerator helix_pairing_ct;
	helix_pairing_ct.element_name( components::HelixPairing::class_name() )
		.complex_type_naming_func ( & pairing_subelement_ct_namer )
		.description( "XRW TO DO" )
		.write_complex_type_to_schema( xsd );


	XMLSchemaComplexTypeGenerator helix_sheet_pairing_ct;
	helix_sheet_pairing_ct.element_name( components::HelixSheetPairing::class_name() )
		.complex_type_naming_func ( & pairing_subelement_ct_namer )
		.description( "XRW TO DO" )
		.write_complex_type_to_schema( xsd );


	XMLSchemaElementOP strand_pairing_subelement( new XMLSchemaElement );
	strand_pairing_subelement->name( components::StrandPairing::class_name() );
	strand_pairing_subelement->type_name( pairing_subelement_ct_namer( components::StrandPairing::class_name() ) );

	XMLSchemaElementOP helix_pairing_subelement( new XMLSchemaElement );
	helix_pairing_subelement->name( components::HelixPairing::class_name() );
	helix_pairing_subelement->type_name( pairing_subelement_ct_namer( components::HelixPairing::class_name() ) );

	XMLSchemaElementOP helix_sheet_pairing_subelement( new XMLSchemaElement );
	helix_sheet_pairing_subelement->name( components::HelixSheetPairing::class_name() );
	helix_sheet_pairing_subelement->type_name( pairing_subelement_ct_namer( components::HelixSheetPairing::class_name() ) );

	XMLSchemaModelGroupOP pairing_choice( new XMLSchemaModelGroup );
	pairing_choice->type( xsmgt_choice );
	pairing_choice->append_particle( helix_pairing_subelement );
	pairing_choice->append_particle( strand_pairing_subelement );
	pairing_choice->append_particle( helix_sheet_pairing_subelement );

	XMLSchemaModelGroup pairing_group;
	pairing_group.group_name( pairing_group_namer() );
	pairing_group.append_particle( pairing_choice );
	xsd.add_top_level_element( pairing_group );
}







void
CompoundArchitect::parse_connection_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	connection::ConnectionArchitectOP architect;
	if ( tag->getName() == "Connection" ) {
		architect = connection::ConnectionArchitectOP(
			new connection::ConnectionArchitect( tag->getOption< std::string >( "name" ) ) );
	} else {
		std::stringstream msg;
		msg << "Bad connection architect type in tag: " << *tag << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
	architect->parse_my_tag( tag, data );
	if ( ! architect->id().empty() ) connection::store_connection_architect( architect, data );
	connections_.push_back( architect );
}

void
CompoundArchitect::parse_pairing_tags( utility::tag::Tag const & tag )
{
	for ( utility::tag::Tag::tags_t::const_iterator t=tag.getTags().begin(); t!=tag.getTags().end(); ++t ) {
		parse_pairing_tag( **t );
	}
}

void
CompoundArchitect::parse_pairing_tag( utility::tag::Tag const & tag )
{
	components::SegmentPairingOP pair = components::SegmentPairing::create( tag.getName() );
	pair->parse_my_tag( tag );
	pairings_.push_back( pair );
}

CompoundArchitect::StructureDataOP
CompoundArchitect::design( core::pose::Pose const & pose, core::Real & random ) const
{
	StructureDataOP sd( new StructureData( id() ) );
	for ( DeNovoArchitectCOPs::const_iterator a=architects_.begin(); a!=architects_.end(); ++a ) {
		StructureDataOP sub_sd = (*a)->design( pose, random );
		if ( !sub_sd ) {
			std::stringstream msg;
			msg << class_name() << "::design(): Architect " << (*a)->id()
				<< " did not produce a StructureData object. Current SD ="
				<< *sd << std::endl;
			utility_exit_with_message( msg.str() );
		}
		if ( sub_sd->segments_begin() == sub_sd->segments_end() ) {
			std::stringstream msg;
			msg << class_name() << "::design(): Architect " << (*a)->id()
				<< " produced a StructureData object without segments. Current SD ="
				<< *sd << " Architect-produced SD = " << *sub_sd << std::endl;
			utility_exit_with_message( msg.str() );
		}
		sd->merge( *sub_sd );
	}

	for ( connection::ConnectionArchitectCOPs::const_iterator c=connections_.begin(); c!=connections_.end(); ++c ) {
		(*c)->apply( *sd, random );
	}

	for ( components::SegmentPairingCOPs::const_iterator p=pairings_.begin(); p!=pairings_.end(); ++p ) {
		sd->add_pairing( **p );
	}
	return sd;
}

void
CompoundArchitect::add_architect( DeNovoArchitect const & architect )
{
	architects_.push_back( architect.clone() );
}

void
CompoundArchitect::add_connection( connection::ConnectionArchitect const & connection )
{
	connections_.push_back( connection.clone() );
}

///////////////////////////////////////////////////////////////////////////////
/// Creator
///////////////////////////////////////////////////////////////////////////////
std::string
CompoundArchitectCreator::keyname() const
{
	return CompoundArchitect::class_name();
}

DeNovoArchitectOP
CompoundArchitectCreator::create_architect( std::string const & architect_id ) const
{
	return DeNovoArchitectOP( new CompoundArchitect( architect_id ) );
}

} //protocols
} //denovo_design
} //architects
