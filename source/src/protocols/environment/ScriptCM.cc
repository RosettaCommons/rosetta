// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/environment/ScriptCM.cc
/// @author Justin R. Porter
/// @author Brian D. Weitzner
/// @author Oliver F. Lange

// Unit Headers
#include <protocols/environment/ScriptCM.hh>
#include <protocols/environment/ScriptCMCreator.hh>

// Package headers

// Project headers
#include <core/environment/DofPassport.hh>
#include <protocols/environment/DofUnlock.hh>

#include <core/kinematics/MoveMap.hh>

#include <protocols/moves/MoverFactory.hh>
#include <protocols/moves/NullMover.hh>

//Utility Headers
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>

// tracer
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#ifdef WIN32
#include <basic/datacache/WriteableCacheableMap.hh>
#endif

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace environment {

std::string const NOT_SET = "[NOT_SET]";

static basic::Tracer tr( "protocols.environment.ScriptCM", basic::t_info );

using namespace core::environment;
using namespace protocols::environment;

// creator
// XRW TEMP std::string
// XRW TEMP ScriptCMCreator::keyname() const {
// XRW TEMP  return ScriptCM::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ScriptCMCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ScriptCM );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ScriptCM::mover_name() {
// XRW TEMP  return "ScriptCM";
// XRW TEMP }

ScriptCM::ScriptCM():
	ClientMover(),
	name_( NOT_SET ),
	client_( new moves::NullMover() )
{}

void ScriptCM::passport_updated(){
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap() );

	if ( has_passport() ) {
		passport()->render_movemap( mm );
	} else {
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_jump( true );
	}

	client()->set_movemap( mm );

}

void ScriptCM::initialize( core::pose::Pose& pose ){
	DofUnlock activation( pose.conformation(), passport() );
	client()->initialize( pose );
}

void ScriptCM::apply( core::pose::Pose& pose ){
	DofUnlock activation( pose.conformation(), passport() );
	client()->apply( pose );
}

void ScriptCM::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const& filters,
	protocols::moves::Movers_map const& mover_map,
	core::pose::Pose const& pose ) {
	name_ = tag->getOption< std::string >( "name" );

	for ( auto tag_it = tag->getTags().begin();
			tag_it != tag->getTags().end(); ++tag_it ) {
		TagCOP subtag = *tag_it;
		if ( (*tag_it)->getName() == "Mover" ) {
			std::string const& client_name = subtag->getOption< std::string >( "name" );
			tr.Debug << " Interpreting tag with name " << subtag->getName() << " as existing mover with name '"
				<< client_name << "'" << std::endl;
			if ( mover_map.find( client_name ) != mover_map.end() ) {
				set_client( mover_map.find( client_name )->second );
			} else {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Undefined mover '"+client_name+"'." );
			}
		} else if ( claims::EnvClaim::is_claim( subtag->getName() ) ) {
			tr.Debug << " Interpreting tag with name " << subtag->getName() << " as new claim." << std::endl;
			add_claim( claims::EnvClaim::make_claim( subtag->getName(), utility::pointer::static_pointer_cast< ClientMover >( get_self_ptr() ), subtag, datamap ) );
		} else {
			tr.Debug << " Interpreting tag with name " << subtag->getName() << " as a new mover." << std::endl;

			set_client( moves::MoverFactory::get_instance()->newMover( subtag, datamap, filters, mover_map, pose ) );

			if ( subtag->hasOption( "name" ) ) {
				tr.Warning << "Mover " << subtag->getOption< std::string >( "name" ) << " will not be availiable to"
					<< " reference by name. It will exist only within " << this->get_name() << std::endl;
			}
		}
	}
}

void ScriptCM::set_client( moves::MoverOP mover_in ) {
	moves::MoveMapMoverOP mover_ptr = utility::pointer::dynamic_pointer_cast< moves::MoveMapMover > ( mover_in );

	if ( !mover_ptr ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "The "+mover_in->type()+" named '"+mover_in->get_name()+
			"' doesn't implement MoveMapMover and can't be used by the ScriptCM." );
	}

	client_ = mover_ptr;
}

void ScriptCM::add_claim( claims::EnvClaimOP claim ) {
	tr.Debug << "  " << this->get_name() << " added EnvClaim " << *claim << "added." << std::endl;
	claim_list_.push_back( claim );
}


claims::EnvClaims ScriptCM::yield_claims( core::pose::Pose const&,
	basic::datacache::WriteableCacheableMapOP ){
	tr.Debug << this->get_name() << " yielding " << claim_list_.size() << " EnvClaims." << std::endl;
	return claim_list_;
}

// XRW TEMP std::string ScriptCM::get_name() const {
// XRW TEMP  return "ScriptCM("+name()+")";
// XRW TEMP }

moves::MoverOP ScriptCM::fresh_instance() const {
	return ClientMoverOP( new ScriptCM() );
}

moves::MoverOP ScriptCM::clone() const{
	return ClientMoverOP( new ScriptCM( *this ) );
}

std::string ScriptCM::get_name() const {
	return mover_name();
}

std::string ScriptCM::mover_name() {
	return "ScriptCM";
}

void ScriptCM::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute(
		"name", xs_string,
		"This mover requires its name attribute to be specified.");


	//So this can actually take three types of subelement
	//1) tag with name "Mover" (already defined)
	//2) An EnvClaim of some sort
	//3) Some other mover

	//So, we actually have to define a group for this

	protocols::moves::MoverFactory::get_instance()->define_mover_xml_schema( xsd );
	protocols::environment::claims::EnvClaim::define_envclaim_schema_group( xsd );

	AttributeList mover_tag_attlist;
	mover_tag_attlist + XMLSchemaAttribute::required_attribute(
		"name", xs_string,
		"Name of a previously defined mover to set as a client mover");


	XMLSchemaComplexTypeGenerator mover_tag_ct;
	mover_tag_ct.element_name( "Mover" )
		.complex_type_naming_func ( & scriptcm_subelement_ct_namer )
		.description( "Specify a previously defined mover to use as a client mover" )
		.add_attributes( mover_tag_attlist )
		.write_complex_type_to_schema( xsd );


	XMLSchemaElementOP mover_tag_subelement( new XMLSchemaElement );
	mover_tag_subelement->name( "Mover" );
	mover_tag_subelement->type_name( scriptcm_subelement_ct_namer( "Mover" ));

	XMLSchemaModelGroupOP mover_reference_group( new XMLSchemaModelGroup(  moves::MoverFactory::mover_xml_schema_group_name() ) );


	XMLSchemaModelGroupOP envclaim_reference_group( new XMLSchemaModelGroup(  claims::EnvClaim::envclaim_group_name() ) );

	XMLSchemaModelGroupOP scriptcm_subelement_choice( new XMLSchemaModelGroup );
	scriptcm_subelement_choice->type( xsmgt_choice );
	scriptcm_subelement_choice->append_particle( mover_reference_group );
	scriptcm_subelement_choice->append_particle( envclaim_reference_group );
	XMLSchemaModelGroup scriptcm_subelement_group;
	scriptcm_subelement_group.group_name( scriptcm_group_name() )
		.append_particle( scriptcm_subelement_choice );
	xsd.add_top_level_element( scriptcm_subelement_group );


	XMLSchemaSimpleSubelementList ssl;
	ssl.add_group_subelement ( & scriptcm_group_name );


	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"XSD_XRW: TO DO",
		attlist, ssl );
}

std::string
ScriptCM::scriptcm_group_name(){
	return "scriptcm";
}

std::string
ScriptCM::scriptcm_subelement_ct_namer( std::string tag_name ){
	return "scriptcm_" + tag_name + "_complex_type";
}



std::string ScriptCMCreator::keyname() const {
	return ScriptCM::mover_name();
}

protocols::moves::MoverOP
ScriptCMCreator::create_mover() const {
	return protocols::moves::MoverOP( new ScriptCM );
}

void ScriptCMCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ScriptCM::provide_xml_schema( xsd );
}


} // environment
} // protocols
