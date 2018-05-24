// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SegmentContactOrderScorer.cc
///
/// @brief Favors assemblies whose segments form contacts with segments distant in the assembly
/// @author Sharon Guffy

//Unit headers
#include <protocols/sewing/scoring/SegmentContactOrderScorer.hh>
#include <protocols/sewing/scoring/SegmentContactOrderScorerCreator.hh>
//Project headers
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
#include <protocols/sewing/scoring/AssemblyScorerFactory.hh>
//Core headers
#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/pose/motif/reference_frames.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/scoring/motif/util.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <numeric/xyzTransform.hh>
#include <utility/io/ozstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
namespace protocols {
namespace sewing  {
namespace scoring {

static basic::Tracer TR("protocols.sewing.scoring.SegmentContactOrderScorer");

SegmentContactOrderScorer::SegmentContactOrderScorer():
	MotifScorer()
{
}

SegmentContactOrderScorer::SegmentContactOrderScorer( SegmentContactOrderScorer const & src ):
	MotifScorer( src ),
	weight_( src.weight_ )
{
}

///@details Use the negative normalized motif score
core::Real
SegmentContactOrderScorer::score(
	data_storage::SmartAssemblyCOP assembly
) {
	this->set_old_last_score( this->get_last_score() );
	core::Real nscore = this->contact_order_score( assembly );
	this->set_last_score( nscore );
	return this->get_last_score();
}

core::Real
SegmentContactOrderScorer::contact_order_score( data_storage::SmartAssemblyCOP assembly ){
	//Iterate over segments
	core::Size n_counter = 0;
	core::Size c_counter = 0;
	data_storage::SmartSegmentOP current_segment = assembly->get_n_terminal_segment();
	core::Real total_score = 0;
	//We'll have a single initial iteration that just finds out how many of the segments in the assembly we're actually counting; this will be what we use to normalize and what we use to compare vs counter in the equation
	core::Size num_helices = 0;
	while ( current_segment != nullptr ) {
		if ( current_segment->get_dssp_code() != 'L' ) {
			++num_helices;
		}
		current_segment = current_segment->get_c_terminal_neighbor();
	}
	current_segment = assembly->get_n_terminal_segment();
	while ( current_segment != nullptr ) {
		//Iterate over the assembly starting at the N terminus until you find a contact and make a note of how many segments you had to go through to get there (if we reach the segment we're looking at, just return that number)
		data_storage::SmartSegmentOP n_contact_segment = assembly->get_n_terminal_segment();
		while ( n_contact_segment != nullptr ) {
			if ( n_contact_segment->get_segment_id() == current_segment->get_segment_id() ) {
				break;
			}
			//Determine whether these segments are in contact. If they are, then break out of the loop.
			bool found_contact = false;
			for ( data_storage::SmartSewingResidueOP current_res : current_segment->get_residue_vector() ) {
				if ( found_contact ) {
					break;
				}
				char ss1 = current_segment->get_dssp_code();
				numeric::xyzTransform< core::Real > stub1 = get_stub( current_res );
				char aa1 = res_type_set_->name_map( current_res->get_amino_acid_type() ).name1();
				for ( data_storage::SmartSewingResidueOP second_res : n_contact_segment->get_residue_vector() ) {
					if ( found_contact ) {
						break;
					}
					char ss2 = n_contact_segment->get_dssp_code();
					numeric::xyzTransform< core::Real > stub2 = get_stub( second_res );
					char aa2 = res_type_set_->name_map( second_res->get_amino_acid_type() ).name1();
					if ( get_score( stub1, ss1, aa1, stub2, ss2, aa2 ) != 0 ) {
						found_contact = true;
					}
				}
			}
			if ( found_contact ) {
				break;
			}
			++n_counter;
			n_contact_segment = n_contact_segment->get_c_terminal_neighbor();
		}
		//Do the same for the C terminus
		data_storage::SmartSegmentOP c_contact_segment = assembly->get_c_terminal_segment();
		while ( c_contact_segment != nullptr ) {
			if ( c_contact_segment->get_segment_id() == current_segment->get_segment_id() ) {
				break;
			}
			//Determine whether these segments are in contact. If they are, then break out of the loop.
			bool found_contact = false;
			for ( data_storage::SmartSewingResidueOP current_res : current_segment->get_residue_vector() ) {
				if ( found_contact ) {
					break;
				}
				char ss1 = current_segment->get_dssp_code();
				numeric::xyzTransform< core::Real > stub1 = get_stub( current_res );
				char aa1 = res_type_set_->name_map( current_res->get_amino_acid_type() ).name1();
				for ( data_storage::SmartSewingResidueOP second_res : c_contact_segment->get_residue_vector() ) {
					if ( found_contact ) {
						break;
					}
					char ss2 = c_contact_segment->get_dssp_code();
					numeric::xyzTransform< core::Real > stub2 = get_stub( second_res );
					char aa2 = res_type_set_->name_map( second_res->get_amino_acid_type() ).name1();
					if ( get_score( stub1, ss1, aa1, stub2, ss2, aa2 ) != 0 ) {
						found_contact = true;
					}
				}
			}
			if ( found_contact ) {
				break;
			}
			++c_counter;
			c_contact_segment = c_contact_segment->get_n_terminal_neighbor();
		}

		//Have some function that converts each of those numbers into a score

		total_score = total_score + ( n_counter * n_counter ) + (c_counter * c_counter );


	}

	//Normalize the score by the number of SSEs calculated earlier
	return total_score / num_helices;
}


void
SegmentContactOrderScorer::set_options_from_tag(
	utility::tag::TagCOP scorer,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const &)
{
	//No data members to set!

	if ( !scorer->hasOption( "weight" ) ) {
		TR << "No weight provided for scorer " << scorer->getName() << " Using default value" << std::endl;
	}

	weight_ = scorer->getOption<core::Real>( "weight", 1.0 );

	TR.Debug << "Created SegmentContactOrderScorer from tag!" << std::endl;

}


core::Real
SegmentContactOrderScorer::get_weight() const{
	return weight_;
}

void
SegmentContactOrderScorer::set_weight( core::Real weight ){
	weight_ = weight;
}

void
SegmentContactOrderScorer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	//only attribute is the weight

	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "How heavily will this term be weighted during scoring?", "1.0"  );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyScorerFactory::assembly_scorer_ct_namer )
		.element_name( SegmentContactOrderScorer::type_name() )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.description( "Favors assemblies whose segments form contacts with segments distant in the assembly" )
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations
}

std::string
SegmentContactOrderScorer::type_name(){
	return "SegmentContactOrderScorer";
}



//////Creator methods//////////
AssemblyScorerOP
SegmentContactOrderScorerCreator::create_assembly_scorer() const{
	return AssemblyScorerOP( new SegmentContactOrderScorer );
}



std::string
SegmentContactOrderScorerCreator::keyname() const{
	return SegmentContactOrderScorer::type_name();
}

void
SegmentContactOrderScorerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	SegmentContactOrderScorer::provide_xml_schema( xsd );


}

/////End Creator methods///////


} //scoring namespace
} //sewing namespace
} //protocols namespace
