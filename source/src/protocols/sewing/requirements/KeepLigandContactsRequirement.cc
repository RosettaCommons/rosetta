// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/KeepLigandContactsRequirement.cc
/// @brief a Requirement that an Assembly have less than a certain number of clashes
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/sewing/requirements/KeepLigandContactsRequirement.hh>
#include <protocols/sewing/requirements/KeepLigandContactsRequirementCreator.hh>
#include <protocols/sewing/requirements/AssemblyRequirementFactory.hh>
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/LigandResidue.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/conformation/Atom.hh>
#include <numeric/xyzVector.hh>
static basic::Tracer TR( "protocols.sewing.requirements.KeepLigandContactsRequirement" );


namespace protocols {
namespace sewing {
namespace requirements {

KeepLigandContactsRequirement::KeepLigandContactsRequirement():
	AssemblyRequirement()
{
	contact_distance_cutoff_ = 2.5;
}

KeepLigandContactsRequirement::~KeepLigandContactsRequirement(){}

KeepLigandContactsRequirement::KeepLigandContactsRequirement( KeepLigandContactsRequirement const & src ):
	AssemblyRequirement( src )
{

	contact_distance_cutoff_ = src.contact_distance_cutoff_;

}



KeepLigandContactsRequirementOP
KeepLigandContactsRequirement::clone() const {
	return KeepLigandContactsRequirementOP( new KeepLigandContactsRequirement( *this ) );
}
std::pair<bool,bool>
KeepLigandContactsRequirement::test(data_storage::SmartAssemblyOP assembly) {
	//Both false if it fails, both true if it passes
	core::Real distance = 0;
	numeric::xyzVector< core::Real > coord1;
	numeric::xyzVector< core::Real > coord2;
	//First iterate over all ligands
	for ( std::pair< core::Size, data_storage::LigandResidueOP > ligand_pair: assembly->get_local_ligands() ) {
		//For each ligand, iterate over contacts
		for ( data_storage::LigandContactOP contact: ligand_pair.second->get_current_contacts() ) {
			if ( contact->segment_id == 0 || contact->residue_number == 0 ) {
				continue;
			}
			//Get the xyz for the ligand atom
			coord1 = ligand_pair.second->get_atom( contact->ligand_atom ).xyz();
			//Get the xyz for the residue atom
			coord2 = assembly->get_segment( contact->segment_id )->get_residue( contact->residue_number )->get_atom( contact->residue_atom ).xyz();
			//Measure the distance between them
			distance = coord1.distance( coord2 );
			//If the distance is greater than the contact distance cutoff, return ( false, false )
			if ( distance > contact_distance_cutoff_ ) {
				return std::make_pair( false, false );
			}
		}
	}
	return std::make_pair( true, true );
}


void
KeepLigandContactsRequirement::set_options_from_tag(
	utility::tag::TagCOP requirement_tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	TR << "Setting up KeepLigandContactsRequirement" << std::endl;
	//maximum_contacts_to_lose_ = requirement_tag->getOption< core::Size >( "maximum_contacts_to_lose", 0 );
	//TR << "Number of lost contacts allowed: " << maximum_contacts_to_lose_ << std::endl;
	contact_distance_cutoff_ = requirement_tag->getOption< core::Real >( "contact_distance_cutoff", 2.5 );
	TR << "Maximum distance between contact atoms (Angstroms): " << contact_distance_cutoff_ << std::endl;
}

void
KeepLigandContactsRequirement::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	AttributeList attributes = KeepLigandContactsRequirement::get_xml_attributes();
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyRequirementFactory::assembly_requirement_ct_namer )
		.element_name( KeepLigandContactsRequirement::type_name() )
		.add_attributes( attributes )
		.description( "Fails if an assembly's ligands lose more than a set number of contacts" )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations
}

utility::tag::AttributeList
KeepLigandContactsRequirement::get_xml_attributes(){
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		//+ XMLSchemaAttribute::attribute_w_default( "maximum_contacts_to_lose", xsct_non_negative_integer, "Maximum number of total ligand contacts that can be lost during assembly before the requirement fails", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "contact_distance_cutoff", xsct_real, "Maximum distance between two contact atoms before the contact is considered broken", "2.5" );
	return attributes;
}

core::Real
KeepLigandContactsRequirement::get_contact_distance_cutoff() const{
	return contact_distance_cutoff_;
}


void
KeepLigandContactsRequirement::set_contact_distance_cutoff( core::Real setting ){
	contact_distance_cutoff_ = setting;
}
//////Creator methods//////////

void
KeepLigandContactsRequirementCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	KeepLigandContactsRequirement::provide_xml_schema( xsd );
}



AssemblyRequirementOP
KeepLigandContactsRequirementCreator::create_requirement() const{
	return AssemblyRequirementOP( new KeepLigandContactsRequirement() );
}

std::string
KeepLigandContactsRequirement::type_name(){
	return "KeepLigandContactsRequirement";
}




std::string
KeepLigandContactsRequirementCreator::keyname() const{
	return KeepLigandContactsRequirement::type_name();
}
/////End Creator methods///////







} //protocols
} //sewing
} //requirements






