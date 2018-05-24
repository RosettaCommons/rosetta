// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/LengthInResiduesRequirement.cc
/// @brief a Requirement that an Assembly be within a certain range of lengths
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/sewing/requirements/LengthInResiduesRequirement.hh>
#include <protocols/sewing/requirements/LengthInResiduesRequirementCreator.hh>
#include <protocols/sewing/requirements/AssemblyRequirementFactory.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
static basic::Tracer TR( "protocols.sewing.requirements.LengthInResiduesRequirement" );


namespace protocols {
namespace sewing {
namespace requirements {

LengthInResiduesRequirement::LengthInResiduesRequirement():
	AssemblyRequirement(),
	minimum_length_( 0 ),
	maximum_length_( 10000 )
{
	test_results_.first = true;
	test_results_.second = true;
}

LengthInResiduesRequirement::LengthInResiduesRequirement(core::Size min_length, core::Size max_length):
	AssemblyRequirement(),
	minimum_length_( min_length ),
	maximum_length_( max_length )
{
	test_results_.first = true;
	test_results_.second = true;
}
LengthInResiduesRequirement::~LengthInResiduesRequirement(){}

LengthInResiduesRequirement::LengthInResiduesRequirement( LengthInResiduesRequirement const & src):
	AssemblyRequirement( src ),
	minimum_length_( src.minimum_length_),
	maximum_length_( src.maximum_length_ )
{
}



LengthInResiduesRequirementOP
LengthInResiduesRequirement::clone() const {
	return LengthInResiduesRequirementOP( new LengthInResiduesRequirement( *this ) );
}

std::pair<bool,bool>
LengthInResiduesRequirement::test(data_storage::SmartAssemblyOP assembly){
	length_ = assembly->get_length();
	test_results_.first = true;
	test_results_.second = false;
	if ( length_ > minimum_length_ ) {
		test_results_.second = true;
		if ( length_ > maximum_length_ ) {
			test_results_.first = false;
			test_results_.second = false;
		}
	}
	return test_results_;
}


//Getters and Setters
core::Size
LengthInResiduesRequirement::get_minimum_length() const{
	return minimum_length_;
}

core::Size
LengthInResiduesRequirement::get_maximum_length() const{
	return maximum_length_;
}

void
LengthInResiduesRequirement::set_minimum_length( core::Size setting ){
	minimum_length_ = setting;
}

void
LengthInResiduesRequirement::set_maximum_length( core::Size setting ){
	maximum_length_ = setting;
}


void
LengthInResiduesRequirement::set_options_from_tag(
	utility::tag::TagCOP requirement_tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	TR << "Setting up LengthInResiduesRequirement" << std::endl;
	maximum_length_ = requirement_tag->getOption< core::Size >( "maximum_length", 10000 );
	TR << "Maximum length: " << maximum_length_ << std::endl;
	minimum_length_ = requirement_tag->getOption< core::Size >( "minimum_length", 0 );
	TR << "Minimum length: " << minimum_length_ << std::endl;
}

void
LengthInResiduesRequirement::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "maximum_length", xsct_non_negative_integer, "Maximum number of residues to allow in the assembly", "10000"  )
		+ XMLSchemaAttribute::attribute_w_default( "minimum_length", xsct_non_negative_integer,  "Minimum number of residues in the final assembly", "0" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyRequirementFactory::assembly_requirement_ct_namer )
		.element_name( LengthInResiduesRequirement::type_name() )
		.add_attributes( attributes )
		.description( "Checks the number of segments in the assembly" )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations

}


//////Creator methods//////////

void
LengthInResiduesRequirementCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{

	LengthInResiduesRequirement::provide_xml_schema( xsd );

}

AssemblyRequirementOP
LengthInResiduesRequirementCreator::create_requirement() const{
	return AssemblyRequirementOP( new LengthInResiduesRequirement() );
}

std::string
LengthInResiduesRequirement::type_name(){
	return "LengthInResiduesRequirement";
}

std::string
LengthInResiduesRequirementCreator::keyname() const{
	return LengthInResiduesRequirement::type_name();
}
/////End Creator methods///////


} //protocols
} //sewing
} //requirements






