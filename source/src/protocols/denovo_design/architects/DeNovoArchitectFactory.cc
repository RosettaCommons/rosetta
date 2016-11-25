// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/DeNovoArchitectFactory.cc
/// @brief Creates DeNovo architects
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.hh>

// Architect creators
#include <protocols/denovo_design/architects/BetaSheetArchitectCreator.hh>
#include <protocols/denovo_design/architects/BlueprintArchitectCreator.hh>
#include <protocols/denovo_design/architects/CompoundArchitectCreator.hh>
#include <protocols/denovo_design/architects/HelixArchitectCreator.hh>
#include <protocols/denovo_design/architects/PoseArchitectCreator.hh>
#include <protocols/denovo_design/architects/StrandArchitectCreator.hh>

// Protocol headers
#include <protocols/denovo_design/architects/DeNovoArchitect.hh>
#include <protocols/denovo_design/architects/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.architects.DeNovoArchitectFactory" );

namespace protocols {
namespace denovo_design {
namespace architects {

DeNovoArchitectFactory::DeNovoArchitectFactory():
	utility::SingletonBase< DeNovoArchitectFactory >(),
	creators_()
{
	// eventually, use registrators, but for now we can just put the creators here
	add_creator( DeNovoArchitectCreatorOP( new BetaSheetArchitectCreator ) );
	add_creator( DeNovoArchitectCreatorOP( new BlueprintArchitectCreator ) );
	add_creator( DeNovoArchitectCreatorOP( new CompoundArchitectCreator ) );
	add_creator( DeNovoArchitectCreatorOP( new HelixArchitectCreator ) );
	add_creator( DeNovoArchitectCreatorOP( new PoseArchitectCreator ) );
	add_creator( DeNovoArchitectCreatorOP( new StrandArchitectCreator ) );
}

DeNovoArchitectFactory::~DeNovoArchitectFactory()
{}

DeNovoArchitectOP
DeNovoArchitectFactory::create_instance(
	std::string const & architect_name,
	std::string const & architect_id ) const
{
	ArchitectCreatorMap::const_iterator arch = creators_.find( architect_name );
	if ( arch == creators_.end() ) {
		std::stringstream msg;
		msg << "DeNovoArchitectFactory::create_instance(): No architect of type "
			<< architect_name << " is registered with the DeNovoArchitectFactory." << std::endl;
		msg << "ID value was " << architect_id << std::endl;
		msg << "You need to register DeNovoArchitects in protocols/denovo_design/architects/DeNovoArchitectFactory.cc" << std::endl;
		msg << "Registered architects: ";
		for ( ArchitectCreatorMap::const_iterator a=creators_.begin(); a!=creators_.end(); ++a ) {
			if ( a != creators_.begin() ) msg << ", ";
			msg << a->first;
		}
		msg << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return arch->second->create_architect( architect_id );
}

DeNovoArchitectOP
DeNovoArchitectFactory::create_from_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data ) const
{
	std::string const keyname = tag->getName();
	std::string const architect_id = tag->getOption< std::string >( "name", "" );
	DeNovoArchitectOP new_arch = create_instance( keyname, architect_id );
	new_arch->parse_my_tag( tag, data );
	if ( !new_arch->id().empty() ) architects::store_denovo_architect( new_arch, data );
	return new_arch;
}

void
DeNovoArchitectFactory::add_creator( DeNovoArchitectCreatorOP creator )
{
	if ( !creator ) {
		std::stringstream msg;
		msg << "DeNovoArchitectFactory::add_creator(): adding a null creator!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	ArchitectCreatorMap::const_iterator arch = creators_.find( creator->keyname() );
	if ( arch != creators_.end() ) {
		std::stringstream msg;
		msg << "DeNovoArchitectFactory:: trying to register a creator named "
			<< creator->keyname() << " twice!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( !creators_.insert( std::make_pair( creator->keyname(), creator ) ).second ) {
		std::stringstream msg;
		msg << "DeNovoArchitectFactory:: failed inserting creator named "
			<< creator->keyname() << " into the map!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
}



void
DeNovoArchitectFactory::define_architect_group( utility::tag::XMLSchemaDefinition & xsd ) const{
	try{
		utility::tag::define_xml_schema_group(
			creators_,
			architect_group_name(),
			& DeNovoArchitectFactory::complex_type_name_for_architect,
			xsd );
	} catch( utility::excn::EXCN_Msg_Exception const & e ) {
		throw utility::excn::EXCN_Msg_Exception( "Could not generate an XML Schema for Architects from DeNovoArchitectFactory; offending class"
			" must call protocols::denovo_design::architects::complex_type_name_for_architect when defining"
			" its XML Schema\n" + e.msg() );
	}
}

std::string
DeNovoArchitectFactory::architect_group_name(){
	return "denovo_architect";
}

std::string
DeNovoArchitectFactory::complex_type_name_for_architect( std::string const & tag_name ){
	return "denovo_architect_" + tag_name + "_complex_type";
}

void
DeNovoArchitectFactory::xsd_architect_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & architect_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes)
{
	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_architect )
		.element_name( architect_type )
		.description( description )
		.add_attributes( attributes )
		.write_complex_type_to_schema( xsd );
}






} //protocols
} //denovo_design
} //architects
