// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/DsspSpecificLengthRequirement.cc
/// @brief a Requirement that the segments of an Assembly with a specific dssp code be within a certain range of lengths
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/sewing/requirements/DsspSpecificLengthRequirement.hh>
#include <protocols/sewing/requirements/DsspSpecificLengthRequirementCreator.hh>
#include <protocols/sewing/requirements/AssemblyRequirementFactory.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
static basic::Tracer TR( "protocols.sewing.requirements.DsspSpecificLengthRequirement" );


namespace protocols {
namespace sewing {
namespace requirements {

DsspSpecificLengthRequirement::DsspSpecificLengthRequirement():
	AssemblyRequirement(),
	dssp_code_( 'X' ),
	minimum_length_( 0 ),
	maximum_length_( 100 )
{
}


DsspSpecificLengthRequirement::DsspSpecificLengthRequirement( DsspSpecificLengthRequirement const & src ):
	AssemblyRequirement( src ),
	dssp_code_( src.dssp_code_ ),
	minimum_length_( src.minimum_length_ ),
	maximum_length_( src.maximum_length_ )
{
}



DsspSpecificLengthRequirementOP
DsspSpecificLengthRequirement::clone() const {
	return DsspSpecificLengthRequirementOP( new DsspSpecificLengthRequirement( *this ) );
}
std::pair<bool,bool>
DsspSpecificLengthRequirement::test(data_storage::SmartAssemblyOP assembly){
	test_results_.first = true;
	test_results_.second = true;
	current_segment_ = assembly->get_n_terminal_segment()->get_c_terminal_neighbor();
	while ( current_segment_ != nullptr && current_segment_->get_c_terminal_neighbor() != nullptr ) {
		if ( current_segment_->get_dssp_code() == dssp_code_ ) {
			if ( current_segment_->get_length() < minimum_length_ || current_segment_->get_length() > maximum_length_ ) {
				test_results_.first = false;
				test_results_.second = false;
			}
		}
		current_segment_ = current_segment_->get_c_terminal_neighbor();
	}
	return test_results_;
}

void
DsspSpecificLengthRequirement::set_options_from_tag(
	utility::tag::TagCOP requirement_tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	TR << "Setting up DsspSpecificLengthRequirement" << std::endl;
	dssp_code_ = requirement_tag->getOption< char >( "dssp_code", 'X' );
	TR << "DSSP code: " << dssp_code_ << std::endl;
	maximum_length_ = requirement_tag->getOption< core::Size >( "maximum_length", 100 );
	TR << "Maximum length: " << maximum_length_ << std::endl;
	minimum_length_ = requirement_tag->getOption< core::Size >( "minimum_length", 0 );
	TR << "Minimum length: " << minimum_length_ << std::endl;
}

void
DsspSpecificLengthRequirement::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	//Define our enum
	//enum dssp_options { H, L, E, X };

	XMLSchemaRestriction dssp_enum;
	dssp_enum.name( "dssp_enum" );
	dssp_enum.base_type( xs_string );
	dssp_enum.add_restriction( xsr_enumeration, "H" );
	dssp_enum.add_restriction( xsr_enumeration, "L" );
	dssp_enum.add_restriction( xsr_enumeration, "E" );
	dssp_enum.add_restriction( xsr_enumeration, "X" );
	xsd.add_top_level_element( dssp_enum );

	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "dssp_code", "dssp_enum", "DSSP code whose length the requirement is restricting", "X"  )
		+ XMLSchemaAttribute::attribute_w_default( "maximum_length", xsct_non_negative_integer,  "Maximum number of residues in a segment with the given secondary structure", "100" )
		+ XMLSchemaAttribute::attribute_w_default( "minimum_length", xsct_non_negative_integer, "Minimum number of residues in a segment with the given secondary structure", "0"  );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyRequirementFactory::assembly_requirement_ct_namer )
		.element_name( DsspSpecificLengthRequirement::type_name() )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.description( "Restricts the number of residues in segments with the specified DSSP" )
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations

}


//Getters and Setters
char
DsspSpecificLengthRequirement::get_dssp_code() const{
	return dssp_code_;
}

core::Size
DsspSpecificLengthRequirement::get_minimum_length() const{
	return minimum_length_;
}

core::Size
DsspSpecificLengthRequirement::get_maximum_length() const{
	return maximum_length_;
}

void
DsspSpecificLengthRequirement::set_dssp_code( char setting){
	dssp_code_ = setting;
}

void
DsspSpecificLengthRequirement::set_minimum_length( core::Size setting ){
	minimum_length_ = setting;
}

void
DsspSpecificLengthRequirement::set_maximum_length( core::Size setting ){
	maximum_length_ = setting;
}



//////Creator methods//////////
void
DsspSpecificLengthRequirementCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	DsspSpecificLengthRequirement::provide_xml_schema( xsd );
}



AssemblyRequirementOP
DsspSpecificLengthRequirementCreator::create_requirement() const{
	return AssemblyRequirementOP( new DsspSpecificLengthRequirement() );
}

std::string
DsspSpecificLengthRequirement::type_name(){
	return "DsspSpecificLengthRequirement";
}

std::string
DsspSpecificLengthRequirementCreator::keyname() const{
	return DsspSpecificLengthRequirement::type_name();
}
/////End Creator methods///////




} //protocols
} //sewing
} //requirements






