// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/NonTerminalStartingSegmentRequirement.cc
/// @brief a Requirement that an Assembly have less than a certain number of clashes
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/sewing/requirements/NonTerminalStartingSegmentRequirement.hh>
#include <protocols/sewing/requirements/NonTerminalStartingSegmentRequirementCreator.hh>
#include <protocols/sewing/requirements/AssemblyRequirementFactory.hh>
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
static basic::Tracer TR( "protocols.sewing.requirements.NonTerminalStartingSegmentRequirement" );


namespace protocols {
namespace sewing {
namespace requirements {

NonTerminalStartingSegmentRequirement::NonTerminalStartingSegmentRequirement():
	AssemblyRequirement()
{
	test_results_.first = true;
	test_results_.second = false;
}

NonTerminalStartingSegmentRequirement::~NonTerminalStartingSegmentRequirement()=default;

NonTerminalStartingSegmentRequirement::NonTerminalStartingSegmentRequirement( NonTerminalStartingSegmentRequirement const & src):
	AssemblyRequirement( src )
{
}



NonTerminalStartingSegmentRequirementOP
NonTerminalStartingSegmentRequirement::clone() const {
	return NonTerminalStartingSegmentRequirementOP( new NonTerminalStartingSegmentRequirement( *this ) );
}
std::pair<bool,bool>
NonTerminalStartingSegmentRequirement::test(data_storage::SmartAssemblyOP assembly) {
	//There will never be a case where this requirement is irrevocably false (building more onto the assembly will never be a problem)
	test_results_.first = true;
	//If any vital segments are found in the assembly, the requirement is not yet satisfied

	//This requirement is really simple--just check that the terminal segments are not vital
	if ( assembly->get_n_terminal_segment()->is_vital() || assembly->get_c_terminal_segment()->is_vital() ) {
		test_results_.second = false;
	} else {
		test_results_.second = true;
	}
	return test_results_;
}


void
NonTerminalStartingSegmentRequirement::set_options_from_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	TR << "Setting up NonTerminalStartingSegmentRequirement" << std::endl;
}

void
NonTerminalStartingSegmentRequirement::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	//There are no attributes!
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyRequirementFactory::assembly_requirement_ct_namer )
		.element_name( NonTerminalStartingSegmentRequirement::type_name() )
		.description( "Requires that any required segments be non-terminal" )
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations
}

//////Creator methods//////////
void
NonTerminalStartingSegmentRequirementCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	NonTerminalStartingSegmentRequirement::provide_xml_schema( xsd );

}




AssemblyRequirementOP
NonTerminalStartingSegmentRequirementCreator::create_requirement() const{
	return AssemblyRequirementOP( new NonTerminalStartingSegmentRequirement() );
}

std::string
NonTerminalStartingSegmentRequirement::type_name(){
	return "NonTerminalStartingSegmentRequirement";
}
std::string
NonTerminalStartingSegmentRequirementCreator::keyname() const{
	return NonTerminalStartingSegmentRequirement::type_name();
}

/////End Creator methods///////







} //protocols
} //sewing
} //requirements






