// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/ClashRequirement.cc
/// @brief a Requirement that an Assembly have less than a certain number of clashes
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/sewing/requirements/ClashRequirement.hh>
#include <protocols/sewing/requirements/ClashRequirementCreator.hh>
#include <protocols/sewing/requirements/AssemblyRequirementFactory.hh>
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/conformation/Atom.hh>
static basic::Tracer TR( "protocols.sewing.requirements.ClashRequirement" );


namespace protocols {
namespace sewing {
namespace requirements {

ClashRequirement::ClashRequirement():
	AssemblyRequirement()
{
	maximum_clashes_allowed_ = 0;
	clash_radius_ = 4.0;
}

ClashRequirement::~ClashRequirement(){}

ClashRequirement::ClashRequirement( ClashRequirement const & src ):
	AssemblyRequirement( src )
{

	maximum_clashes_allowed_ = src.maximum_clashes_allowed_;
	clash_radius_ = src.clash_radius_;

}



ClashRequirementOP
ClashRequirement::clone() const {
	return ClashRequirementOP( new ClashRequirement( *this ) );
}
std::pair<bool,bool>
ClashRequirement::test(data_storage::SmartAssemblyOP assembly) {
	// TR << "Beginning clash check" << std::endl;
	current_clashes_ = 0;
	// if(assembly->get_last_change() == 'D' || assembly->get_last_change() == 'R'){
	//  return test_results_;
	// }
	active_segment_ = assembly->get_n_terminal_segment();
	while ( active_segment_->is_c_terminus_fixed() ) {
		runtime_assert( active_segment_ != nullptr );
		reference_segment_ = active_segment_->get_c_terminal_neighbor();
		while ( reference_segment_->is_c_terminus_fixed() ) {
			reference_segment_ = reference_segment_->get_c_terminal_neighbor(); //First check is actually chimaera (active segment) against its C-terminal neighbor (reference segment)
			for ( active_resnum_ = 1; active_resnum_ <= active_segment_->get_length(); ++active_resnum_ ) {
				active_residue_ = active_segment_->get_residue(active_resnum_);
				if ( active_residue_->get_atom_vector().size() == 0 ) {
					TR << "Residue " << active_resnum_ << " in segment " << active_segment_->get_segment_id() << "reports no atoms!" << std::endl;
				}
				for ( reference_resnum_ = 1; reference_resnum_ <= reference_segment_->get_length(); ++reference_resnum_ ) {
					//NOTE: This check should not be necessary since we aren't checking adjacent segments
					//Do not check adjacent residues at borders between segments!
					//Reference segment is c teriminal neighbor of active segment
					if ( reference_resnum_ == 1 && active_resnum_ == active_segment_->get_length() && active_segment_->get_c_terminal_neighbor()->get_segment_id() == reference_segment_->get_segment_id() ) {
						continue;
					}
					//The reverse should never happen since reference segment is always c-terminal to active segment


					reference_residue_ = reference_segment_->get_residue(reference_resnum_);
					if ( reference_residue_->get_atom_vector().size() == 0 ) {
						TR << "Residue " << reference_resnum_ << " in segment " << reference_segment_->get_segment_id() << "reports no atoms!" << std::endl;
						continue;
					}
					// for( active_atom_number_ = 1; active_atom_number_ <= active_residue_->get_atom_vector().size(); ++active_atom_number_ ){
					// for( reference_atom_number_ = 1; reference_atom_number_ <= reference_residue_->get_atom_vector().size(); ++reference_atom_number_ ){
					//Only checking for CA clashes to avoid false positives
					if ( active_residue_->get_atom( 2 ).xyz().distance(reference_residue_->get_atom( 2 ).xyz()) <= clash_radius_ ) {
						TR << "Self clash detected!" << std::endl;
						current_clashes_++;
						if ( current_clashes_ > maximum_clashes_allowed_ ) {
							//TR << "ClashRequirement failed!" << std::endl;
							test_results_.first = false;
							test_results_.second = false;
							return test_results_;
						}
					}
					// }
					// }
				}
			}
		}


		if ( assembly->get_partner()!=nullptr ) {
			for ( active_resnum_ = 1; active_resnum_ <= active_segment_->get_length(); ++active_resnum_ ) {
				active_residue_ = active_segment_->get_residue(active_resnum_);
				if ( active_residue_->get_atom_vector().size() == 0 ) {
					TR << "Residue " << active_resnum_ << " in segment " << active_segment_->get_segment_id() << "reports no atoms!" << std::endl;
				}
				for ( partner_resnum_ = 1; partner_resnum_ < assembly->get_partner()->total_residue(); partner_resnum_++ ) {
					core::conformation::Atom partner_atom = assembly->get_partner()->residue(partner_resnum_).atoms().at(2);
					if ( active_residue_->get_atom(2).xyz().distance(partner_atom.xyz()) <= clash_radius_ ) {
						//TR << "Partner clash detected!" << std::endl;
						current_clashes_++;
						if ( current_clashes_ > maximum_clashes_allowed_ ) {
							//TR << "ClashRequirement failed!" << std::endl;
							test_results_.first = false;
							test_results_.second = false;
							return test_results_;
						}
					}
				}
			}
		}
		active_segment_ = active_segment_->get_c_terminal_neighbor();
	}
	if ( assembly->get_partner()!=nullptr ) {
		for ( active_resnum_ = 1; active_resnum_ <= active_segment_->get_length(); ++active_resnum_ ) {
			active_residue_ = active_segment_->get_residue(active_resnum_);
			if ( active_residue_->get_atom_vector().size() == 0 ) {
				TR << "Residue " << active_resnum_ << " in segment " << active_segment_->get_segment_id() << "reports no atoms!" << std::endl;
			}
			for ( partner_resnum_ = 1; partner_resnum_ < assembly->get_partner()->total_residue(); partner_resnum_++ ) {
				core::conformation::Atom partner_atom = assembly->get_partner()->residue(partner_resnum_).atoms().at(2);
				if ( active_residue_->get_atom(2).xyz().distance(partner_atom.xyz()) <= clash_radius_ ) {
					//TR << "Partner clash detected!" << std::endl;
					current_clashes_++;
					if ( current_clashes_ > maximum_clashes_allowed_ ) {
						//TR << "ClashRequirement failed!" << std::endl;
						test_results_.first = false;
						test_results_.second = false;
						return test_results_;
					}
				}
			}
		}
	}

	if ( current_clashes_ > maximum_clashes_allowed_ ) {
		//TR << "ClashRequirement failed!" << std::endl;
		test_results_.first = false;
		test_results_.second = false;
	} else {
		//TR << "ClashRequirement passed!" << std::endl;
		test_results_.first = true;
		test_results_.second = true;
	}
	return test_results_;
}


void
ClashRequirement::set_options_from_tag(
	utility::tag::TagCOP requirement_tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	TR << "Setting up ClashRequirement" << std::endl;
	maximum_clashes_allowed_ = requirement_tag->getOption< core::Size >( "maximum_clashes_allowed", 0 );
	TR << "Number of clashes allowed: " << maximum_clashes_allowed_ << std::endl;
	clash_radius_ = requirement_tag->getOption< core::Real >( "clash_radius", 4.0 );
	TR << "Clash radius (Angstroms): " << clash_radius_ << std::endl;
}

void
ClashRequirement::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	/*
	AttributeList attributes;
	attributes
	+ XMLSchemaAttribute::attribute_w_default( "maximum_clashes_allowed", xsct_non_negative_integer, "Maximum number of clashes to allow in the assembly", "0" )
	+ XMLSchemaAttribute::attribute_w_default( "clash_radius", xs_decimal, "Radius in Angstroms within which two residues are considered to be clashing", "4.0" );
	*/
	AttributeList attributes = ClashRequirement::get_xml_attributes();
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyRequirementFactory::assembly_requirement_ct_namer )
		.element_name( ClashRequirement::type_name() )
		.add_attributes( attributes )
		.description( "Checks for clashes between segments in the assembly" )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations
}

utility::tag::AttributeList
ClashRequirement::get_xml_attributes(){
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "maximum_clashes_allowed", xsct_non_negative_integer, "Maximum number of clashes to allow in the assembly", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "clash_radius", xsct_real, "Radius in Angstroms within which two residues are considered to be clashing", "4.0" );
	return attributes;
}

core::Size
ClashRequirement::get_maximum_clashes_allowed() const{
	return maximum_clashes_allowed_;
}

core::Real
ClashRequirement::get_clash_radius() const{
	return clash_radius_;
}


void
ClashRequirement::set_maximum_clashes_allowed( core::Size setting){
	maximum_clashes_allowed_ = setting;
}

void
ClashRequirement::set_clash_radius( core::Real setting ){
	clash_radius_ = setting;
}
//////Creator methods//////////

void
ClashRequirementCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	ClashRequirement::provide_xml_schema( xsd );
}



AssemblyRequirementOP
ClashRequirementCreator::create_requirement() const{
	return AssemblyRequirementOP( new ClashRequirement() );
}

std::string
ClashRequirement::type_name(){
	return "ClashRequirement";
}




std::string
ClashRequirementCreator::keyname() const{
	return ClashRequirement::type_name();
}
/////End Creator methods///////







} //protocols
} //sewing
} //requirements






