// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/SizeInSegmentsRequirement.cc
/// @brief a Requirement that an Assembly be within a certain range of sizes
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/sewing/requirements/SizeInSegmentsRequirement.hh>
#include <protocols/sewing/requirements/SizeInSegmentsRequirementCreator.hh>
#include <protocols/sewing/requirements/AssemblyRequirementFactory.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
static basic::Tracer TR( "protocols.sewing.requirements.SizeInSegmentsRequirement" );


namespace protocols {
namespace sewing {
namespace requirements {

SizeInSegmentsRequirement::SizeInSegmentsRequirement():
	AssemblyRequirement(),
	minimum_size_( 0 ),
	maximum_size_( 10 )
{
	test_results_.first = true;
	test_results_.second = true;
}

SizeInSegmentsRequirement::SizeInSegmentsRequirement(core::Size min_size, core::Size max_size):
	AssemblyRequirement(),
	minimum_size_( min_size ),
	maximum_size_( max_size )
{
	test_results_.first = true;
	test_results_.second = true;
}
SizeInSegmentsRequirement::~SizeInSegmentsRequirement()=default;

SizeInSegmentsRequirement::SizeInSegmentsRequirement( SizeInSegmentsRequirement const & src):
	AssemblyRequirement( src ),
	minimum_size_( src.minimum_size_ ),
	maximum_size_( src.maximum_size_ )
{
	test_results_.first = true;
	test_results_.second = true;
}



SizeInSegmentsRequirementOP
SizeInSegmentsRequirement::clone() const {
	return SizeInSegmentsRequirementOP( new SizeInSegmentsRequirement( *this ) );
}

std::pair<bool,bool>
SizeInSegmentsRequirement::test(data_storage::SmartAssemblyOP assembly){
	size_ = assembly->get_size();
	test_results_.first = true;
	test_results_.second = false;
	if ( size_ > minimum_size_ ) {
		test_results_.second = true;
		if ( size_ > maximum_size_ ) {
			test_results_.first = false;
			test_results_.second = false;
		}
	}
	return test_results_;
}

//Getters and Setters
core::Size
SizeInSegmentsRequirement::get_minimum_size() const{
	return minimum_size_;
}

core::Size
SizeInSegmentsRequirement::get_maximum_size() const{
	return maximum_size_;
}

void
SizeInSegmentsRequirement::set_minimum_size( core::Size setting ){
	minimum_size_ = setting;
}

void
SizeInSegmentsRequirement::set_maximum_size( core::Size setting ){
	maximum_size_ = setting;
}


void
SizeInSegmentsRequirement::set_options_from_tag(
	utility::tag::TagCOP requirement_tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	TR << "Setting up SizeInSegmentsRequirement" << std::endl;
	maximum_size_ = requirement_tag->getOption< core::Size >( "maximum_size", 10000 );
	TR << "Maximum size: " << maximum_size_ << std::endl;
	minimum_size_ = requirement_tag->getOption< core::Size >( "minimum_size", 0 );
	TR << "Minimum size: " << minimum_size_ << std::endl;
}

void
SizeInSegmentsRequirement::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "maximum_size", xsct_non_negative_integer, "Maximum number of secondary structure elements (including loops) to allow in the assembly", "10000"  )
		+ XMLSchemaAttribute::attribute_w_default( "minimum_size", xsct_non_negative_integer,  "Minimum number of secondary structure elements (including loops) in the final assembly", "0" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyRequirementFactory::assembly_requirement_ct_namer )
		.element_name( SizeInSegmentsRequirement::type_name() )
		.add_attributes( attributes )
		.description( "Checks the number of segments in the assembly" )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations

}


//////Creator methods//////////

void
SizeInSegmentsRequirementCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{

	SizeInSegmentsRequirement::provide_xml_schema( xsd );

}

AssemblyRequirementOP
SizeInSegmentsRequirementCreator::create_requirement() const{
	return AssemblyRequirementOP( new SizeInSegmentsRequirement() );
}

std::string
SizeInSegmentsRequirement::type_name(){
	return "SizeInSegmentsRequirement";
}

std::string
SizeInSegmentsRequirementCreator::keyname() const{
	return SizeInSegmentsRequirement::type_name();
}
/////End Creator methods///////


} //protocols
} //sewing
} //requirements






