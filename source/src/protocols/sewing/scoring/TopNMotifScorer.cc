// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TopNMotifScorer.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/sewing/scoring/TopNMotifScorer.hh>
#include <protocols/sewing/scoring/TopNMotifScorerCreator.hh>
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

#include <utility/vector1.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <numeric/xyzTransform.hh>
#include <utility/io/ozstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
namespace protocols {
namespace sewing  {
namespace scoring {

static basic::Tracer TR("protocols.sewing.scoring.TopNMotifScorer");


//Why in the world does this not derive from MotifScorer? Feels like code duplication to me.


TopNMotifScorer::TopNMotifScorer():
	MotifScorer()
{
	scores_.clear();
}

TopNMotifScorer::TopNMotifScorer( TopNMotifScorer const & src ):
	MotifScorer( src ),
	weight_( src.weight_ ),
	scores_( src.scores_ ),
	scores_to_keep_( src.scores_to_keep_ )
{
}

///@details Use the negative normalized motif score
core::Real
TopNMotifScorer::score(
	data_storage::SmartAssemblyCOP assembly
) {
	this->set_old_last_score( this->get_last_score() );
	core::Real nscore = this->norm_motif_score( assembly );
	this->set_last_score(  -1.0 * nscore );
	TR.Debug << "Final motif score: " << this->get_last_score() << std::endl;
	return this->get_last_score();
}


///@details Motif score of the entire Assembly divided by total residue
core::Real
TopNMotifScorer::norm_motif_score(
	data_storage::SmartAssemblyCOP assembly
) {
	core::Real fscore = this->full_motif_score( assembly );
	if ( assembly->get_size() <= 1 ) {
		return 0;
	}
	return fscore / core::Real (assembly->get_size()*(assembly->get_size()-1) );
}


///@details use Will's Motif score to calculate the motif score for the entire Assembly
core::Real
TopNMotifScorer::full_motif_score(
	data_storage::SmartAssemblyCOP assembly
) {
	TR << "starting scoring" << std::endl;
	core::Real score = 0.0;
	core::Real running_score = 0.0;
	utility::vector1<core::Real>::iterator it;
	data_storage::SmartSegmentOP current_segment = assembly->get_n_terminal_segment();
	runtime_assert( current_segment != nullptr );
	//Note that this segment is pulled from the assembly's set of segments, not the master set
	while ( current_segment != nullptr ) {
		//TR << "Scoring segment " << current_segment->get_segment_id() << " of length " << current_segment->get_length() << std::endl;
		//DSSP actually lives in segment and won't change, so just get it here
		char ss1 = current_segment->get_dssp_code();
		//Iterate through the residues in this segment
		//I need the indices to avoid double counting
		//  for( data_storage::SmartSewingResidueOP current_res : current_segment->get_residue_vector() ){
		for ( core::Size current_res_index = 1; current_res_index <=current_segment->get_length() ; ++current_res_index ) {
			//New version of get_stub just takes the SmartSewingResidueOP
			data_storage::SmartSewingResidueCOP current_res = current_segment->get_const_residue_vector().at( current_res_index );
			numeric::xyzTransform< core::Real > stub1 = get_stub( current_res );
			//We already have the DSSP, just get the aa type
			char aa1 = res_type_set_->name_map( current_res->get_amino_acid_type() ).name1();

			//Now loop through all of the segments starting at this segment
			data_storage::SmartSegmentOP second_segment = current_segment;
			while ( second_segment != nullptr ) {
				scores_.clear();
				scores_.resize(scores_to_keep_,0.0);
				//TR << "Scoring against segment " << second_segment->get_segment_id() << " of length " << second_segment->get_length() << std::endl;
				char ss2 = second_segment->get_dssp_code();
				//for( data_storage::SmartSewingResidueOP second_res : second_segment->get_residue_vector() ){
				for ( core::Size second_res_index = 1; second_res_index <= second_segment->get_length(); ++second_res_index ) {
					//Double check that we're not comparing the residue to itself
					if ( current_segment->get_segment_id() == second_segment->get_segment_id() && second_res_index == current_res_index ) { //These two pointers actually point to the same residue
						continue;
					}
					data_storage::SmartSewingResidueCOP second_res = second_segment->get_const_residue_vector().at( second_res_index );
					numeric::xyzTransform< core::Real > stub2 = get_stub( second_res );
					char aa2 = res_type_set_->name_map( second_res->get_amino_acid_type() ).name1();
					//TR << "Scoring atom pair with ss1 " << ss1 << " ss2 " << ss2 << " aa1 " << aa1 << " aa2 " << aa2 << std::endl;
					//TopN mods start here!
					score = get_score(stub1, ss1, aa1, stub2, ss2, aa2);
					if ( score > 0 ) {
						scores_.push_back(score);
					}
					//TopN mods end here!

				} //end for second_res
				second_segment = second_segment->get_c_terminal_neighbor();
				for ( core::Size i = 1; i <= scores_to_keep_; i++ ) {
					it = max_element(std::begin(scores_),std::end(scores_));
					running_score += *it;
					scores_.erase(it);
				}
			} //end while second_segment
			TR << "running score is: " << running_score << std::endl;
		} //end for current_res
		//if( !( current_segment->is_c_terminus_fixed() ) ){
		// break;
		//}
		current_segment = current_segment->get_c_terminal_neighbor();
	}
	return running_score;
}




void
TopNMotifScorer::set_options_from_tag(
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

	if ( scorer->hasOption("scores_to_keep") ) {
		scores_to_keep_ = scorer->getOption<core::Size> ("scores_to_keep", 1);
	}
	TR.Debug << "Created TopNMotifScorer from tag!" << std::endl;

}


std::string
TopNMotifScorer::get_name() const {
	return "TopNMotifScorer";
}

core::Real
TopNMotifScorer::get_weight() const{
	return weight_;
}

void
TopNMotifScorer::set_weight( core::Real weight ){
	weight_ = weight;
}
core::Size
TopNMotifScorer::get_scores_to_keep() const{
	return scores_to_keep_;
}
void
TopNMotifScorer::set_scores_to_keep( core::Size new_scores){
	scores_to_keep_ = new_scores;
}

void
TopNMotifScorer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	//only attribute is the weight

	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "How heavily will this term be weighted during scoring?", "1.0"  )
		+ XMLSchemaAttribute::attribute_w_default( "scores_to_keep", xsct_non_negative_integer, "How many scores from each pair should be counted?", "1");

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyScorerFactory::assembly_scorer_ct_namer )
		.element_name( TopNMotifScorer::type_name() )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.description( "Basic Motif score among all helices" )
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations



}

std::string
TopNMotifScorer::type_name(){
	return "TopNMotifScorer";
}



//////Creator methods//////////
AssemblyScorerOP
TopNMotifScorerCreator::create_assembly_scorer() const{
	return AssemblyScorerOP( new TopNMotifScorer );
}



std::string
TopNMotifScorerCreator::keyname() const{
	return TopNMotifScorer::type_name();
}

void
TopNMotifScorerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	TopNMotifScorer::provide_xml_schema( xsd );


}

/////End Creator methods///////



} //scoring namespace
} //sewing namespace
} //protocols namespace
