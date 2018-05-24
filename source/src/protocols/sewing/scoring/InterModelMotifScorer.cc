// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file InterModelMotifScorer.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/sewing/scoring/InterModelMotifScorer.hh>
#include <protocols/sewing/scoring/InterModelMotifScorerCreator.hh>

//Package headers
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

static basic::Tracer TR("protocols.sewing.scoring.InterModelMotifScorer");

InterModelMotifScorer::InterModelMotifScorer():
	MotifScorer(),
	weight_( 10.0 )
{
}
InterModelMotifScorer::InterModelMotifScorer( InterModelMotifScorer const & src):
	MotifScorer( src ),
	weight_( src.weight_ )
{
}

///@details Use the negative normalized motif score
core::Real
InterModelMotifScorer::score(
	data_storage::SmartAssemblyCOP assembly
) {
	this->set_old_last_score( this->get_last_score() );
	this->set_last_score( -1.0 * full_motif_score(assembly) );
	TR.Debug << "Final intermodel motif score: " << this->get_last_score() << std::endl;
	return this->get_last_score();
}

///@details use Will's Motif score to calculate the motif score for interactions between
///a given segment and segments from other models. Divide by total number of segments
core::Real
InterModelMotifScorer::full_motif_score(
	data_storage::SmartAssemblyCOP assembly
) {
	TR.Debug << "Scoring InterModelMotifScorer"<< std::endl;
	core::Real score = 0.0;
	core::Size counter = 0; //This keeps track of how many interactions we have scored
	data_storage::SmartSegmentOP current_segment = assembly->get_n_terminal_segment();
	//Note that this segment is pulled from the assembly's set of segments, not the master set
	while ( current_segment != nullptr ) {
		//TR << "Current segment ID " << current_segment->get_segment_id() << std::endl;
		//TR << "Current segment dssp code " << current_segment->get_dssp_code() << std::endl;
		char ss1 = current_segment->get_dssp_code();
		//Iterate through the residues in this segment
		for ( data_storage::SmartSewingResidueOP current_res : current_segment->get_residue_vector() ) {
			//TR <<"Scoring a residue" << std::endl;
			//New version of get_stub just takes the SmartSewingResidueOP
			numeric::xyzTransform< core::Real > stub1 = get_stub( current_res );
			//We already have the DSSP, just get the aa type
			char aa1 = res_type_set_->name_map( current_res->get_amino_acid_type() ).name1();
			//We want to score all these interactions once and only once
			//We no longer really need to worry about the residue index though!
			data_storage::SmartSegmentOP second_segment = current_segment;
			//In this case, we definitely don't want to score the segment against itself, so we can go ahead and move down the line
			second_segment = current_segment->get_c_terminal_neighbor();
			while ( second_segment != nullptr ) {
				//TR << "Scoring residue against the next segment " << std::endl;
				//We'll need to check whether second_segment is part of the same model as current_segment
				//As luck would have it, we only need to look backwards (toward the N terminus), not forwards
				//possible ways to accomplish this:

				//1) Start by getting the n-terminal segment for this substructure (not crossing chimarae) and move forward until w
				//get to second_segment. If we don't find current_segment, we're good.
				data_storage::SmartSegmentOP tester = data_storage::SmartSegment::get_n_most_segment(second_segment,false );
				//This is only slightly hacky
				bool same_substructure = false;
				while ( tester->get_segment_id() != second_segment->get_segment_id() && !same_substructure ) {
					//TR << "Testing if the second segment is part of this segment's model" << std::endl;
					if ( tester->get_segment_id() == current_segment->get_segment_id() ) {
						//we don't want to compare current_segment and second_segment
						same_substructure = true;
						break;
					}
					if ( !tester->is_c_terminus_fixed() ) {
						TR << "Something went wrong with InterModelMotifScore's tester segment! Aborting!" << std::endl;
						utility_exit_with_message( "Tester segment ran off the end of the Assembly before finding second_segment!" );
					}
					tester = tester->get_c_terminal_neighbor();
				}
				if ( same_substructure ) {
					//TR <<"Same substructure detected" << std::endl;
					second_segment = second_segment->get_c_terminal_neighbor();
					continue;
				}
				//Okay, in theory, this means that second_segment and current_segment
				//are on different substructures
				//Now we can find the motif score between them like normal
				char ss2 = second_segment->get_dssp_code();
				for ( data_storage::SmartSewingResidueOP second_res : second_segment->get_residue_vector() ) {
					//TR << "Comparing two residues" << std::endl;
					numeric::xyzTransform< core::Real > stub2 = get_stub( second_res );
					char aa2 = res_type_set_->name_map( second_res->get_amino_acid_type() ).name1();
					score += get_score(stub1, ss1, aa1, stub2, ss2, aa2);
					++counter;
				}//for second_res
				second_segment = second_segment->get_c_terminal_neighbor();
			}//while second_segment
		} //for current_res
		//if( !( current_segment->is_c_terminus_fixed() ) ){
		// break;
		//}
		current_segment = current_segment->get_c_terminal_neighbor();
	} //while current_segment
	if ( counter == 0 ) { return score; }
	return score / assembly->get_length();
}


void
InterModelMotifScorer::set_options_from_tag(
	utility::tag::TagCOP scorer,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

	if ( !scorer->hasOption( "weight" ) ) {
		TR << "No weight provided for scorer " << scorer->getName() << " Using default value" << std::endl;
	}

	weight_ = scorer->getOption<core::Real>( "weight", 1.0 );

	//No data members to set!
	TR.Debug << "Created InterModelMotifScorer from tag." << std::endl;

}
std::string
InterModelMotifScorer::get_name() const {
	return "InterModelMotifScorer";
}

core::Real
InterModelMotifScorer::get_weight() const{
	return weight_;
}

void
InterModelMotifScorer::set_weight( core::Real weight ){
	weight_ = weight;
}

void
InterModelMotifScorer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){

	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real,  "How heavily will this term be weighted during scoring?", "1.0" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyScorerFactory::assembly_scorer_ct_namer )
		.element_name( type_name() )
		.description( "Basic Motif score among non-adjacent helices" )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations
}





//////Creator methods//////////
void
InterModelMotifScorerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	InterModelMotifScorer::provide_xml_schema( xsd );

}

AssemblyScorerOP
InterModelMotifScorerCreator::create_assembly_scorer() const{
	return AssemblyScorerOP( new InterModelMotifScorer() );
}

std::string
InterModelMotifScorer::type_name(){
	return "InterModelMotifScorer";
}


std::string
InterModelMotifScorerCreator::keyname() const{
	return InterModelMotifScorer::type_name();
}
/////End Creator methods///////






} //scoring namespace
} //sewing namespace
} //protocols namespace
